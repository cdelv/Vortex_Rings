#include "navier_solver.hpp"
#include <fstream>
#include <string>

using namespace mfem;
using namespace navier;

//Configuration Functions
struct Config
{
    //Numerical Method Parameters
    int n = 6;
    int serial_refinements = 0;
    int parallel_refinements = 1;
    int order = 2;

    //Time Parameters
    int vis_freq = 50;
    double dt = 0.0001;
    double t_final = 0.002;

    //Box Parameters
    double Lx = 2;
    double Ly = 1;
    double Lz = 1;

    //Piston Parameters
    double Piston_R =0.15; //Radius
    double Piston_T=0.05; //Thickness
    double Piston_W=2;    //Omega
    double Piston_P=2*M_PI/Piston_W; //Period

    //Tube Parameters 
    double Tube_R = 0.2; //Radius
    double Tube_L = 0.5; //Lenght
    double Tube_T = 0.05; //Thickness (REVISAR)

    //Piston Outflow Parameters
    double End_R = 0.2;

    //Physical Parameters
    double kinvis = 0.0005;
    double atm_pressure = 0.0;
    double gravity = 9.8;
    double Da = 5e-6;               //Darcy Number, for Solid Obj, 1e-6 < Da < 1e-5
    double l = Tube_L;             // critical lenght 
    double eta = kinvis/(Da*l*l); //-----> kinvis/Da*l^2  Brinkman Penalization

    double Delta = 1/Da;

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
    virtual void Eval(mfem::Vector &a, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
    void SetVel(mfem::GridFunction* gfvel){vel=gfvel;};
    void SetTime(double tt);
    double Piston(const Vector &X, double t);
    double Tube(const Vector &X, double t);
    double End(const Vector &X, double t);
    double Chi(const Vector &X, double t);
    double Impermeability(const Vector &X, double t);
    void Create_Chi_Coefficient(ParGridFunction &CChi);

private:
    mfem::GridFunction* vel = nullptr;
    double t=0.0,x=0.0,vx=0.0;
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
    NavierSolver *flowsolver = nullptr;
    BrinkPenalAccel *piston = nullptr;

    ParMesh *pmesh = new ParMesh();
    {
    //Load Mesh (Pointer to Delete it After Parallel Mesh is Created)
    Mesh mesh = Mesh::MakeCartesian3D(2*Parameters.n, Parameters.n, Parameters.n, Element::QUADRILATERAL, Parameters.Lx, Parameters.Ly, Parameters.Lz);;
    mesh.EnsureNodes();
    int dim = mesh.Dimension();

    //Refine Serial Mesh
    for (int i = 0; i < Parameters.serial_refinements; ++i)
        mesh.UniformRefinement();

    //Make Parallel Mesh
    pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
    }
    
    //Refine Parallel Mesh
    for (int ii = 0; ii < Parameters.parallel_refinements; ii++)
        pmesh->UniformRefinement();

    //Create H1 Element Collections
    H1_FECollection vfec = H1_FECollection(Parameters.order, pmesh->Dimension());
    H1_FECollection pfec = H1_FECollection(Parameters.order);

    //Create Finite Element Spaces
    ParFiniteElementSpace vfes = ParFiniteElementSpace(pmesh, &vfec, pmesh->Dimension()); //Vector Space
    ParFiniteElementSpace pfes = ParFiniteElementSpace(pmesh, &pfec);                   //Scalar Space

    //Define Navier Solver
    flowsolver = new NavierSolver(pmesh, Parameters.order, Parameters.kinvis);
    flowsolver->EnablePA(true);
    flowsolver->EnableNI(true);
    flowsolver->EnableVerbose(false);

    //Tengo algo confundido con las fronteras y la cond inicial
    //Define Velocity Initial Conditions
    ParGridFunction *u_0 = flowsolver->GetCurrentVelocity();
    VectorFunctionCoefficient u_bdr(pmesh->Dimension(), Initial_Velocity);
    u_0->ProjectCoefficient(u_bdr);

    //Define Pressure Initial Conditions
    ParGridFunction *p_0 = flowsolver->GetCurrentPressure();
    ConstantCoefficient  p_bdr(Parameters.atm_pressure);
    p_0->ProjectCoefficient(p_bdr);

    //Define Velocity Boundary Conditions
    Array<int> attr(pmesh->bdr_attributes.Max());
    //Inflow      Outflow        Side down      Side up       Side left      Side right
    attr[4] = 1;  attr[2] = 0;   attr[0] = 0;   attr[5] = 0;  attr[1] = 0;   attr[3] = 0;
    flowsolver->AddVelDirichletBC(Vel_Boundary_Condition, attr);

    //Define  Pressure Boundary Conditions
    Array<int> attr2(pmesh->bdr_attributes.Max());
    //Inflow       Outflow         Side down       Side up        Side left       Side right
    attr2[4] = 0;  attr2[2] = 1;   attr2[0] = 1;   attr2[5] = 1;  attr2[1] = 1;   attr2[3] = 1;
    flowsolver->AddPresDirichletBC(Press_Boundary_Condition, attr2);


    //Define Solution Pointers 
    ParGridFunction *u = flowsolver->GetCurrentVelocity(); //Velocity Solution
    ParGridFunction *p = flowsolver->GetCurrentPressure(); //Pressure Solution

    //Calculate Vorticity
    ParGridFunction w = ParGridFunction(&vfes);
    CurlGridFunctionCoefficient curl_u(u);
    w.ProjectCoefficient(curl_u);

    //Create Domain Attr Array
    Array<int> domain_attr(pmesh->bdr_attributes.Max());
    domain_attr=1;

    //Add Gravity Acceleration Term
    //flowsolver.AddAccelTerm(Gravity,domain_attr);  

    //Create Solid object
    piston = new BrinkPenalAccel(pmesh->Dimension());
    piston->SetVel(flowsolver->GetCurrentVelocity());
    piston->SetTime(t);
    flowsolver->AddAccelTerm(piston,domain_attr);   

    //Calculate Object Contour
    ParGridFunction Chi = ParGridFunction(&pfes);
    piston->Create_Chi_Coefficient(Chi);

    //Set Up Solver (Must be After Adding Accel Terms)
    flowsolver->Setup(Parameters.dt);

    //Paraview Visualization
    ParaViewDataCollection paraview_out = ParaViewDataCollection("results/graph", pmesh);
    paraview_out.SetLevelsOfDetail(Parameters.order);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.RegisterField("Object", &Chi);
    paraview_out.RegisterField("Pressure", p);
    paraview_out.RegisterField("Velocity", u);
    paraview_out.RegisterField("Vorticity", &w);

    if(mpi.Root())
        std::cout << "step" << "\t" << "t" << "\t" << "dt" << "\t" << "print" << "\n";

    //Perform Time Integration
    for (int step = 0; !last_step; ++step)
    {
        //Check if Last Step
        if (t + Parameters.dt >= Parameters.t_final - Parameters.dt / 2)
            last_step = true;

        //Time Step
        flowsolver->Step(t, Parameters.dt, step);

        //Update Time of The Velocity Boundary Condition
        u_bdr.SetTime(t);

        //Update The Solid Object
        piston->SetVel(flowsolver->GetCurrentVelocity());
        piston->SetTime(t);

        if(mpi.Root()){
            std::cout.flush();
            std::cout << step << "\t" << t << "\t" << Parameters.dt << "\t" << vis_print << "\r";
        }

        //Print Data for Visualization
        if (step%Parameters.vis_freq==0)
        {   
            CurlGridFunctionCoefficient u_curl(u);
            w.ProjectCoefficient(u_curl);
            piston->Create_Chi_Coefficient(Chi);
            paraview_out.SetCycle(vis_print);
            paraview_out.SetTime(t);
            paraview_out.Save();
            vis_print++;
        }
    }

    flowsolver->PrintTimingData();

    //with MPI sesion theres no neet to Finalize MPI
    //Free memory, Dont have any pointer yet

    delete pmesh;
    delete flowsolver;

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
    double eta = Parameters.eta;// Impermeability(X,t);

    //Get Fluid Velocity
    vel->GetVectorValue(T,ip,U);

    //Only de Piston Should Move
    U0(0)=vx*Piston(X,t);

    //The - Sing is Already Included
    a(0)=chi*eta*(U0(0)-U(0));
    a(1)=chi*eta*(U0(1)-U(1));
    a(2)=chi*eta*(U0(2)-U(2));
}
void BrinkPenalAccel::SetTime(double tt)
{
	//Update Time
	t=tt;

	//Update Piston Position and Velocity
	double dist = 1.5;
	if (t<=Parameters.Piston_P/4)
	{
	    x=Parameters.Tube_L*(dist+std::sin(Parameters.Piston_W*t+3*M_PI_2));
	    vx=Parameters.Tube_L*Parameters.Piston_W*std::cos(Parameters.Piston_W*t+3*M_PI_2);
	}
	else{
		x=Parameters.Tube_L*dist;
		vx=0;
	}
}
double BrinkPenalAccel::Piston(const Vector &X, double t)
{   
    double r = std::pow(X(1)-Parameters.Ly/2,2)+std::pow(X(2)-Parameters.Lz/2,2);
    double c = -0.5*std::tanh(Parameters.Delta*(r-std::pow(Parameters.Piston_R,2)))+0.5; //Disk Radius
    double h = -0.5*std::tanh(Parameters.Delta*(std::abs(X(0)-x)-Parameters.Piston_T))+0.5; //Disk Thickness
    return c*h;
}
double BrinkPenalAccel::Tube(const Vector &X, double t)
{   
	double mean = (Parameters.Piston_R+Parameters.Tube_R)*0.5;
    double r = std::sqrt(std::pow(X(1)-Parameters.Ly/2,2)+std::pow(X(2)-Parameters.Lz/2,2));
    double c = -0.5*std::tanh(Parameters.Delta*(std::abs(r-mean)-(Parameters.Tube_R-Parameters.Piston_R)*0.5))+0.5; //Wall thickness
    double h = -0.5*std::tanh(Parameters.Delta*(X(0)-Parameters.Tube_L))+0.5; //Tube Lenght
    return c*h;
}
double BrinkPenalAccel::End(const Vector &X, double t)
{   
	double mean = (Parameters.End_R+Parameters.Tube_R)*0.5;
    double r = std::sqrt(std::pow(X(1)-Parameters.Ly/2,2)+std::pow(X(2)-Parameters.Lz/2,2)); 
    double c = -0.5*std::tanh(Parameters.Delta*(std::abs(r-mean)-(Parameters.Tube_R-Parameters.End_R)*0.5))+0.5; //End Opnening
    double h = -0.5*std::tanh(Parameters.Delta*(std::abs(X(0)-Parameters.Tube_L)-Parameters.Piston_T))+0.5; //Face thickness
    return c*h;
}
double BrinkPenalAccel::Chi(const Vector &X, double t)
{
    return Piston(X,t);//+Tube(X,t)+End(X,t);
}
//Inverse of the brinkman penalization permeability
double BrinkPenalAccel::Impermeability(const Vector &X, double t){
    double Phase = 1-Chi(X,t);
    //return Epsilon + pow(Chi(X,t), 2)/(pow(1-Chi(X,t), 3) + Epsilon);
    return Parameters.Da + pow(1-Phase, 2)/(pow(Phase, 3) + Parameters.Da);
} 

void BrinkPenalAccel::Create_Chi_Coefficient(ParGridFunction &CChi)
{   
    //FunctionCoefficient Require a Static Function
    auto Create_Chi = [=](const Vector &X, double t) mutable {return Chi(X,t);};
    FunctionCoefficient chi(Create_Chi);
    CChi.ProjectCoefficient(chi);
}
