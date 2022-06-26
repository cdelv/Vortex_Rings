#include "navier_solver.hpp"
#include <fstream>
#include <string>

using namespace mfem;
using namespace navier;

//Configuration Functions
struct Config
{
    //Numerical Method Parameters
    int serial_refinements = 1;
    int parallel_refinements = 1;
    int order = 2;

    //Time Parameters
    int vis_freq = 50;
    double dt = 0.0001;
    double t_final = 3;

    //Box Parameters
    double Lx = 2;
    double Ly = 1;
    double Lz = 1;

    //Piston Parameters
    double Piston_R =0.15;
    double Piston_T=0.05;
    double Piston_W=2;
    double Piston_P=2*M_PI/Piston_W;

    //Tube Parameters
    double Tube_R = 0.2;
    double Tube_L = 0.5;
    double Tube_T = 0.05;

    //Piston Outflow Parameters
    double End_R = 0.2;

    //Physical Parameters
    double kinvis = 0.0005;
    double atm_pressure = 0.0;
    double gravity = 9.8;
    double Da = 5e-6;               //Darcy Number, for Solid Obj, 1e-6 < Da < 1e-5
    double l = Tube_L;             // critical lenght 
    double eta = kinvis/(Da*l*l); //-----> kinvis/Da*l^2  Brinkman Penalization

    //Read Parameters
    void init(const std::string &parameters_file);
} Parameters;
double get_next_parameter(std::ifstream &file);

//Assembly functions
void Initial_Velocity(const Vector &x, double t, Vector &u);
void Initial_Vorticity(const Vector &x, double t, Vector &u);
void Vel_Boundary_Condition(const Vector &x, double t, Vector &u);
double Press_Boundary_Condition(const Vector &x, double t);
void Gravity(const Vector &x, double t, Vector &a);

//Brinkman term class
class BrinkPenalAccel:public mfem::VectorCoefficient{
public:
    BrinkPenalAccel(int dim):mfem::VectorCoefficient(dim){};
    virtual ~BrinkPenalAccel(){};
    void SetUp(double X,double Y,double Z,double VX,double VY,double VZ,double R,double H);
    void SetVel(mfem::GridFunction* gfvel){vel=gfvel;};
    void SetTime(double tt){t=tt;};
    virtual void Eval(mfem::Vector &a, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
    double Chi(const Vector &X, double t);
    void Create_Chi_Coefficient(ParGridFunction &CChi);
    bool Piston(const Vector &X, double t);

private:
    mfem::GridFunction* vel = nullptr;
    double t=0.0,x=0.0;
};

//Main function
int main(int argc, char *argv[])
{   
    //Init MPI
    MPI_Session mpi(argc, argv);

    //Read Parameters From File
    //std::string parameters_file = "settings/parameters.txt";
    //Parameters.init(parameters_file);

    //Aux Variables
    double t = 0.0;
    int vis_print = 0;
    bool last_step = false;

    //Load Mesh (Pointer to Delete it After Parallel Mesh is Created)
    Mesh *mesh = new Mesh("mesh.msh");
    mesh->EnsureNodes();
    int dim = mesh->Dimension();

    //Refine Serial Mesh
    for (int i = 0; i < Parameters.serial_refinements; ++i)
        mesh->UniformRefinement();

    //Make Parallel Mesh
    ParMesh pmesh = ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    
    //Refine Parallel Mesh
    for (int ii = 0; ii < Parameters.parallel_refinements; ii++)
        pmesh.UniformRefinement();

    //Create H1 Element Collections
    H1_FECollection vfec = H1_FECollection(Parameters.order, pmesh.Dimension());
    H1_FECollection pfec = H1_FECollection(Parameters.order);

    //Create Finite Element Spaces
    ParFiniteElementSpace vfes = ParFiniteElementSpace(&pmesh, &vfec, pmesh.Dimension()); //Vector Space
    ParFiniteElementSpace pfes = ParFiniteElementSpace(&pmesh, &pfec);                   //Scalar Space

    //Calculate Vorticity
    ParGridFunction w = ParGridFunction(&vfes);
    ParGridFunction u = ParGridFunction(&vfes);
    ParGridFunction psi = ParGridFunction(&vfes);
    VectorFunctionCoefficient w_bdr(pmesh.Dimension(), Initial_Vorticity);
    w.ProjectCoefficient(w_bdr);

    {
        ParBilinearForm a(&vfes);
        a.AddDomainIntegrator(new VectorCurlCurlIntegrator);
        a.Assemble();

        ParLinearForm b(&vfes);
        b.AddDomainIntegrator(new VectorDomainLFIntegrator(w_bdr));
        b.Assemble();

        psi = 0.;

        Array<int> attr_psi(pmesh.bdr_attributes.Max());
        attr_psi = 0;   attr_psi[0] = 1;    //Inflow    
        Array<int> tdof_psi;
        vfes.GetEssentialTrueDofs(attr_psi, tdof_psi);

        HypreParMatrix A;
        HypreParVector B(&vfes);
        HypreParVector Psi(&vfes);
        a.FormLinearSystem(tdof_psi, psi, b, A, Psi, B);
        
        HypreBoomerAMG *amg = new HypreBoomerAMG(A);
        amg->SetPrintLevel(0);
    
        //Solve the linear system Ax=B
        HyprePCG *pcg = new HyprePCG(A);
        pcg->SetPreconditioner(*amg);
        pcg->SetPrintLevel(0);
        pcg->SetTol(1e-12);
        pcg->SetMaxIter(200);
        pcg->Mult(B, Psi);
    
        //Recover the solution on each proccesor
        a.RecoverFEMSolution(Psi, b, psi);
        CurlGridFunctionCoefficient curl_psi(&psi);
        u.ProjectCoefficient(curl_psi);
    
        //Delete used memory
        delete amg;
        delete pcg;
    }

    //Paraview Visualization
    ParaViewDataCollection paraview_out = ParaViewDataCollection("results/graph", &pmesh);
    paraview_out.SetLevelsOfDetail(Parameters.order);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.RegisterField("Vorticity", &w);
    paraview_out.RegisterField("Psi", &psi);
    paraview_out.RegisterField("Velocity", &u);
    paraview_out.SetCycle(vis_print);
    paraview_out.SetTime(t);
    paraview_out.Save();

    return 0;
}
//Configuration functions
/*void Config::init(const std::string &parameters_file){
    std::ifstream fparams(parameters_file);
    serial_refinements = (int) get_next_parameter(fparams);
    parallel_refinements = (int) get_next_parameter(fparams);
    order = (int) get_next_parameter(fparams);
    vis_freq = (int) get_next_parameter(fparams);
    dt = get_next_parameter(fparams);
    t_final = get_next_parameter(fparams);

    R_ring = get_next_parameter(fparams);
    ring_center_y = get_next_parameter(fparams);
    ring_center_z = get_next_parameter(fparams);
    V0_ring = get_next_parameter(fparams);
    I_vel_exponent = get_next_parameter(fparams);
    I_active_dts = get_next_parameter(fparams);

    kinvis = get_next_parameter(fparams);
    atm_pressure = get_next_parameter(fparams);
    gravity = get_next_parameter(fparams);
    fparams.close();
}
double get_next_parameter(std::ifstream &file){
    std::string name_parameter;
    getline(file,name_parameter,' ');
    std::string string_value;
    getline(file,string_value);
    return stod(string_value);
}*/
//Assemble functions
void Initial_Velocity(const Vector &x, double t, Vector &u)
{
   double xi = x(0);
   double yi = x(1);
   double zi = x(2);

    u(0) = 0.0;
    u(1) = 0.0;
    u(2) = 0.0;
}

void LinealVortex(const double t, Vector &u){
    double R = 0.3;
    double y0 = 0.5;
    double z0 = 0.5;
    u(0) = 1.;
    u(1) = y0 + R*std::cos(t);
    u(2) = z0 + R*std::sin(t);
}

void TangentVortex(const double t, Vector &u){
    u(0) = 0.;
    u(1) = -std::sin(t);
    u(2) = std::cos(t);
}

void Initial_Vorticity(const Vector &x, double t, Vector &u)
{

    double a = 0.1;
    double y0 = 0.5;
    double z0 = 0.5;

    double xi = x(0);
    double yi = x(1);
    double zi = x(2);
    double theta = std::atan2(zi-z0, yi-y0);
   
    LinealVortex(theta, u);
    double r = u.DistanceTo(x);

    if (r <= a){
        TangentVortex(theta, u);
        u *= 1-r/a;
    } else
        u = 0.;
}
void Vel_Boundary_Condition(const Vector &x, double t, Vector &u)
{
    double xi = x(0);
    double yi = x(1);
    double zi = x(2);

    u(0) = 0.0;
    u(1) = 0.0;
    u(2) = 0.0;
}
double Press_Boundary_Condition(const Vector &x, double t)
{
    return Parameters.atm_pressure;
}
void Gravity(const Vector &x, double t, Vector &a)
{
   double xi = x(0);
   double yi = x(1);
   double zi = x(2);

    a(0) = 0.0;
    a(1) = 0.0;
    a(2) = -Parameters.gravity;
}
//Brinkman term class
void BrinkPenalAccel::Eval(mfem::Vector &a, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
    //Body velocity    Fluid velocity     Fluid coordinates
    Vector U0;         Vector U;          Vector X;
    
    //Set Body and Fluid Velocity and Acceleration 
    U0.SetSize(GetVDim()); U0(0)=0; U0(1)=0; U0(2)=0;
    U.SetSize(GetVDim());  U(0)=0;  U(1)=0;  U(2)=0;
    a.SetSize(GetVDim());  a(0)=0;  a(1)=0;  a(2)=0;

    //Get the physical coordinates of integration point
    T.Transform(ip,X); 

    //Check For the object
    double chi = Chi(X,t);

    //Get Fluid Velocity
    vel->GetVectorValue(T,ip,U);

    double dist = 1.3;
    if (Piston(X,t) && t<=Parameters.Piston_P/4)
    {
        x=Parameters.Tube_L*(dist+std::sin(Parameters.Piston_W*t+3*M_PI_2));
        U0(0)=Parameters.Tube_L*Parameters.Piston_W*std::cos(Parameters.Piston_W*t+3*M_PI_2);
    }
    else
    	x=Parameters.Tube_L*dist;

    //The - Sing is Already Included
    a(0)=chi*Parameters.eta*(U0(0)-U(0));
    a(1)=chi*Parameters.eta*(U0(1)-U(1));
    a(2)=chi*Parameters.eta*(U0(2)-U(2));
}
double BrinkPenalAccel::Chi(const Vector &X, double t)
{
    if(Piston(X,t))
        return 1;

    return 0;
}
void BrinkPenalAccel::Create_Chi_Coefficient(ParGridFunction &CChi)
{   
    //FunctionCoefficient Require a Static Function
    auto Create_Chi = [=](const Vector &X, double t) mutable {return Chi(X,t);};
    FunctionCoefficient chi(Create_Chi);
    CChi.ProjectCoefficient(chi);
}
bool BrinkPenalAccel::Piston(const Vector &X, double t)
{   
    if(std::pow(X(1)-Parameters.Ly/2,2)+std::pow(X(2)-Parameters.Lz/2,2)<std::pow(Parameters.Piston_R,2) && std::abs(X(0)-x) < Parameters.Piston_T)
        return true;

    return false;
}
