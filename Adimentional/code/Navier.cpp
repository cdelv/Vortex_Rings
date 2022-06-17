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

    //Time Parameters (seconds)
    int vis_freq = 50;
    double dt = 0.0001;
    double t_final = 20.00;

    //Piston Parameters (meters)
    double r_int = 0.025;
    double r_ext = 0.05;
    double lenght = 0.5;
    double thickness = 0.05;

    //Piston stroke (rad/s)
    double omega = 100;

    //Piston opening (circular) (meters)
    double r_open = 0.02; 

    //Box parameters (meters)
    double Lx=1.0;
    double Ly=0.5;
    double Lz=0.5;

    //Physical Parameters
    double kinvis = 0.00001488;  //m^2 s^-1 (air viscosity)
    double atm_pressure = 0.0;	// Pa -- m^-1 s^-2 kg
    double gravity = 9.8;	   //m s^-2
    double Da = 5e-6;         //Darcy Number, for Solid Obj, 1e-6 < Da < 1e-5.

    //(fluid density is assumed as 1)
    double ct,cl,eta;

    //Read Parameters
    void init(const std::string &parameters_file);

    //Remove Units
    void Adimentionalize(void);
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
    void SetVel(mfem::GridFunction* gfvel){vel=gfvel;};
    void SetTime(double tt){t=tt;};
    virtual void Eval(mfem::Vector &a, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
    double Chi(const Vector &X, double t);
    void Create_Chi_Coefficient(ParGridFunction &CChi);

private:
    mfem::GridFunction* vel = nullptr;
    double t=0, pos_x=0,pos_aux=0,vel_aux=1;
    bool Piston_Wall(const Vector &X, double t);
    bool Piston(const Vector &X, double t);
    bool Piston_Exit(const Vector &X, double t);
};

//Main function
int main(int argc, char *argv[])
{   
    //Init MPI
    MPI_Session mpi(argc, argv);

    //Read Parameters From File
    //std::string parameters_file = "settings/parameters.txt";
    //Parameters.init(parameters_file);

    //Remove Units
    Parameters.Adimentionalize();

    //Parameters
    double t = 0.0;
    int vis_print = 0;
    bool last_step = false;

    //Load Mesh (Pointer to Delete it After Parallel Mesh is Created)
    //Mesh *mesh = new Mesh("mesh.msh");
    //mesh->EnsureNodes();
    //int dim = mesh->Dimension();

    //                                     nx,             ny,           nz,           elements,               Lx,            Ly,            Lz
    Mesh mesh = Mesh::MakeCartesian3D(2*Parameters.n, Parameters.n, Parameters.n, Element::QUADRILATERAL, Parameters.Lx, Parameters.Ly, Parameters.Lz);
    mesh.EnsureNodes();
    int dim = mesh.Dimension();

    //Refine Serial Mesh
    for (int i = 0; i < Parameters.serial_refinements; ++i)
        mesh.UniformRefinement();

    //Make Parallel Mesh
    ParMesh pmesh = ParMesh(MPI_COMM_WORLD, mesh);
    //delete mesh;
    
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
    /*std::ifstream fparams(parameters_file);
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
    fparams.close();*/
}
double get_next_parameter(std::ifstream &file){
    std::string name_parameter;
    getline(file,name_parameter,' ');
    std::string string_value;
    getline(file,string_value);
    return stod(string_value);
}
void Config::Adimentionalize(void)
{
	//Critical Lenght Will be max(Ly,Lz) for the Analog With Poiseuille Flow
	cl = r_int;

	//Critical Time is the Piston Stroke Time
	ct = 2*M_PI/(4*omega);

	//Time Parameters (seconds)
    vis_freq /= ct;
    dt /= ct;
    t_final /= ct;

    //Piston Parameters (meters)
    r_int /= cl;
    r_ext /= cl;
    lenght /= cl;
    thickness /= cl;

    //Piston stroke (rad/s)
    omega /= ct;

    //Piston opening (circular) (meters)
    r_open /= cl;

    //Box parameters (meters)
    Lx /= cl;
    Ly /= cl;
    Lz /= cl;

    //Physical Parameters
    kinvis /= (cl*cl/ct);         //m^2 s^-1 (air viscosity)
    atm_pressure /= 1/(cl*ct*ct);// Pa --> m^-1 s^-2 kg
    gravity /= (cl/(ct*ct));	   //m s^-2
    eta = kinvis/(Da*cl*cl);   //Kinvis/Da*l^2  Brinkman Penalization.
}
//Assemble functions
void Initial_Velocity(const Vector &x, double t, Vector &u)
{
    u(0) = 0.0;
    u(1) = 0.0;
    u(2) = 0.0;
}
void Vel_Boundary_Condition(const Vector &x, double t, Vector &u)
{
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
    U0.SetSize(GetVDim()); U0(1)=0; U0(2)=0;
    U.SetSize(GetVDim());
    a.SetSize(GetVDim());

    //Armonic motion to the right
    pos_x = Parameters.lenght + Parameters.lenght*std::sin(Parameters.omega*t+3*M_PI/2) + Parameters.thickness/2+pos_aux;
    U0(0) = Parameters.lenght*Parameters.omega*std::cos(Parameters.omega*t+3*M_PI/2)*vel_aux;
    if(pos_x >= Parameters.lenght - Parameters.thickness/2)
    {	
        U0(0)=0;
    	pos_x = Parameters.lenght - Parameters.thickness/2;
        pos_aux = 2*(Parameters.lenght - Parameters.thickness/2);
        vel_aux = 0;
    }

    //Get the physical coordinates at integration point
    T.Transform(ip,X); 

    //Get the fluid velocity at the integration point
    vel->GetVectorValue(T,ip,U);

    //Calculate if the object is present at the integration point
    double chi = Chi(X,t);

    //If piston wall or exit, v_object=0
    if (Piston_Wall(X,t) || Piston_Exit(X,t))
    	U0(0)=0;

    //The - Sing is Already Included
    a(0)=Parameters.eta*(U0(0)-U(0))*chi;
    a(1)=Parameters.eta*(U0(1)-U(1))*chi;
    a(2)=Parameters.eta*(U0(2)-U(2))*chi;
}
double BrinkPenalAccel::Chi(const Vector &X, double t)
{	
    if (Piston_Wall(X,t) || Piston(X,t) || Piston_Exit(X,t))
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
bool BrinkPenalAccel::Piston_Wall(const Vector &X, double t)
{
	if (std::pow(X(1)-Parameters.Ly/2,2)+std::pow(X(2)-Parameters.Lz/2,2)<Parameters.r_ext && std::pow(X(1)-Parameters.Ly/2,2)+std::pow(X(2)-Parameters.Lz/2,2)>Parameters.r_int && X(0)<Parameters.lenght)
		return true;

	return false;
}

bool BrinkPenalAccel::Piston(const Vector &X, double t)
{
	if (std::pow(X(1)-Parameters.Ly/2,2)+std::pow(X(2)-Parameters.Lz/2,2)<Parameters.r_int && X(0)< pos_x + Parameters.thickness/2 )
		return true;

	return false;
}
bool BrinkPenalAccel::Piston_Exit(const Vector &X, double t)
{
	if (std::pow(X(1)-Parameters.Ly/2,2)+std::pow(X(2)-Parameters.Lz/2,2)<Parameters.r_ext && std::pow(X(1)-Parameters.Ly/2,2)+std::pow(X(2)-Parameters.Lz/2,2)>Parameters.r_open && X(0)>Parameters.lenght-Parameters.thickness/2 && X(0)<Parameters.lenght+Parameters.thickness/2)
		return true;

	return false;
}