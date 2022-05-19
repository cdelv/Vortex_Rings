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
    double t_final = 20.00;

    //Initial Conditions Parameters
    double V0_ring = 0.0;
    double R_ring = 0.2;
    double ring_center_y = 0.5;
    double ring_center_z = 0.5;
    double I_vel_exponent = 50.0;
    double I_active_dts = 1.0;

    //Physical Parameters
    double kinvis = 0.0005;
    double atm_pressure = 0.0;
    double gravity = 9.8;
    double Da = 5e-6;               //Darcy Number, for Solid Obj, 1e-6 < Da < 1e-5
    double l = 2;                  // critical lenght (lenght of mesh)
    double eta = kinvis/(Da*l*l); //-----> kinvis/Da*l^2  Brinkman Penalization

    //Read Parameters
    void init(const std::string &parameters_file);
} Parameters;
double get_next_parameter(std::ifstream &file);

//Assembly functions
void Initial_Velocity(const Vector &x, double t, Vector &u);
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

private:
    mfem::GridFunction* vel = nullptr;
    double t=0,x=0.6,y=0.5,z=0.5,vx=0.1,vy=0,vz=0,r=0.4,h=0.2;
    double x0=0, y0=0, z0=0;
};

//Main function
int main(int argc, char *argv[])
{   
    //Init MPI
    MPI_Session mpi(argc, argv);

    //Read Parameters From File
    //std::string parameters_file = "settings/parameters.txt";
    //Parameters.init(parameters_file);

    //Parameters
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

    //Define Navier Solver
    NavierSolver flowsolver(&pmesh, Parameters.order, Parameters.kinvis);
    flowsolver.EnablePA(true);
    flowsolver.EnableNI(true);
    flowsolver.EnableVerbose(true);

    //Define Velocity Initial Conditions
    ParGridFunction *u_0 = flowsolver.GetCurrentVelocity();
    VectorFunctionCoefficient u_bdr(pmesh.Dimension(), Initial_Velocity);
    u_0->ProjectCoefficient(u_bdr);

    //Define Pressure Initial Conditions
    ParGridFunction *p_0 = flowsolver.GetCurrentPressure();
    ConstantCoefficient  p_bdr(Parameters.atm_pressure);
    p_0->ProjectCoefficient(p_bdr);

    //Define Velocity Boundary Conditions
    Array<int> attr(pmesh.bdr_attributes.Max());
    //Inflow       Sides           Outflow
    attr[0] = 1;   attr[1] = 0;    attr[2] = 0;
    flowsolver.AddVelDirichletBC(Vel_Boundary_Condition, attr);

    //Define  Pressure Boundary Conditions
    Array<int> attr2(pmesh.bdr_attributes.Max());
    //Inflow       Sides           Outflow
    attr2[0] = 0;  attr2[1] = 1;   attr2[2] = 1;
    flowsolver.AddPresDirichletBC(Press_Boundary_Condition, attr2);

    //Define Solution Pointers 
    ParGridFunction *u = flowsolver.GetCurrentVelocity(); //Velocity Solution
    ParGridFunction *p = flowsolver.GetCurrentPressure(); //Pressure Solution

    //Calculate Vorticity
    ParGridFunction w = ParGridFunction(&vfes);
    CurlGridFunctionCoefficient curl_u(u);
    w.ProjectCoefficient(curl_u);

    //Create Domain Attr Array
    Array<int> domain_attr(pmesh.bdr_attributes.Max()); 
    domain_attr=1;

    //Add Gravity Acceleration Term
    //flowsolver.AddAccelTerm(Gravity,domain_attr);  

    //Create Solid object
    BrinkPenalAccel piston(pmesh.Dimension());
    piston.SetUp(0.3,0.5,0.5,0.5,0,0,0.2,0.05);
    piston.SetVel(flowsolver.GetCurrentVelocity());
    piston.SetTime(t);
    flowsolver.AddAccelTerm(&piston,domain_attr);   

    //Calculate Object Contour
    ParGridFunction Chi = ParGridFunction(&pfes);
    piston.Create_Chi_Coefficient(Chi);

    //Set Up Solver (Must be After Adding Accel Terms)
    flowsolver.Setup(Parameters.dt);

    //Paraview Visualization
    ParaViewDataCollection paraview_out = ParaViewDataCollection("results/graph", &pmesh);
    paraview_out.SetLevelsOfDetail(Parameters.order);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.RegisterField("Object", &Chi);
    paraview_out.RegisterField("Pressure", p);
    paraview_out.RegisterField("Velocity", u);
    paraview_out.RegisterField("Vorticity", &w);
    paraview_out.SetCycle(vis_print);
    paraview_out.SetTime(t);
    paraview_out.Save();

    //Perform Time Integration
    for (int step = 0; !last_step; ++step)
    {
        //Check if Last Step
        if (t + Parameters.dt >= Parameters.t_final - Parameters.dt / 2)
            last_step = true;

        //Time Step
        flowsolver.Step(t, Parameters.dt, step);

        //Update Time of The Velocity Boundary Condition
        u_bdr.SetTime(t);

        //Update The Solid Object
        piston.SetVel(flowsolver.GetCurrentVelocity());
        piston.SetTime(t);

        //Compute CFL Condition, Must be <= 1
        /*double cfl = flowsolver.ComputeCFL(*u, Parameters.dt);

        if(mpi.Root()){
            std::cout << "step" << "\t" << "t" << "\t" << "dt" << "\t" << "cfl" << "\n";
            std::cout << step << "\t" << t << "\t" << Parameters.dt << "\t" << cfl << "\n";
        }*/

        //Print Data for Visualization
        if (step%Parameters.vis_freq==0)
        {   
            vis_print++;
            CurlGridFunctionCoefficient u_curl(u);
            w.ProjectCoefficient(u_curl);
            piston.Create_Chi_Coefficient(Chi);
            paraview_out.SetCycle(vis_print);
            paraview_out.SetTime(t);
            paraview_out.Save();
        }
    }

    flowsolver.PrintTimingData();

    //with MPI sesion theres no neet to Finalize MPI
    //Free memory, Dont have any pointer yet

    return 0;
}
//Configuration functions
void Config::init(const std::string &parameters_file){
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
}
//Assemble functions
void Initial_Velocity(const Vector &x, double t, Vector &u)
{
   double xi = x(0);
   double yi = x(1);
   double zi = x(2);

    u(0) = 0.0;
    u(1) = 0.0;
    u(2) = 0.0;

  if(pow(yi-Parameters.ring_center_y,2)+pow(zi-Parameters.ring_center_z,2)<pow(Parameters.R_ring,2)){
        u(0) = Parameters.V0_ring*exp(-Parameters.I_vel_exponent*xi);
  }
}
void Vel_Boundary_Condition(const Vector &x, double t, Vector &u)
{
    double xi = x(0);
    double yi = x(1);
    double zi = x(2);

    u(0) = 0.0;
    u(1) = 0.0;
    u(2) = 0.0;

    if(pow(yi-Parameters.ring_center_y,2)+pow(zi-Parameters.ring_center_z,2)<pow(Parameters.R_ring,2) && t<Parameters.I_active_dts*Parameters.dt){
        u(0) = Parameters.V0_ring;
    }
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
void BrinkPenalAccel::SetUp(double X,double Y,double Z,double VX,double VY,double VZ,double R,double H)
{
    x0=X;  y0=Y;  z0=Z;
    x=X;   y=Y;   z=Z;
    vx=VX; vy=VY; vz=VZ;
    r=R;   h=H;
}
void BrinkPenalAccel::Eval(mfem::Vector &a, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
    //Body velocity    Fluid velocity     Fluid coordinates
    Vector U0;         Vector U;          Vector X;
    
    //Set Body and Fluid Velocity and Acceleration 
    U0.SetSize(GetVDim()); U0(0)=vx; U0(1)=vy; U0(2)=vz;
    U.SetSize(GetVDim());  U(0)=0;  U(1)=0;   U(2)=0;
    a.SetSize(GetVDim());  a(0)=0;  a(1)=0;   a(2)=0;

    //Get the physical coordinates of integration point
    T.Transform(ip,X); 

    x=x0+vx*std::tanh(10*t)*t;

    vel->GetVectorValue(T,ip,U);
    double chi = Chi(X,t);

    //The - Sing is Already Included
    a(0)=Parameters.eta*(U0(0)-U(0)*std::tanh(10*t))*chi;
    a(1)=Parameters.eta*(U0(1)-U(1)*std::tanh(10*t))*chi;
    a(2)=Parameters.eta*(U0(2)-U(2)*std::tanh(10*t))*chi;
}
double BrinkPenalAccel::Chi(const Vector &X, double t)
{
    if(std::pow(X(1)-y,2)+std::pow(X(2)-z,2)<std::pow(r,2) && std::abs(X(0)-x) < h)
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