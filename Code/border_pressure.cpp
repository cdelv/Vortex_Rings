#include "../Navier/navier_solver.hpp"
#include <fstream>
#include <string>

using namespace mfem;
using namespace navier;

struct Navier_Parameters
{
    //Parameters
    int serial_refinements = 1;
    int parallel_refinements = 0;
    int order = 2;
    int vis_freq = 50;
    double dt = 0.001;
    double t_final = 0.005;

    double R_ring = 0.2;
    double ring_center_y = 0.5;
    double ring_center_z = 0.5;
    double V0_ring = 7.0;
    double I_vel_exponent = 30.0;
    double I_active_dts = 1.0;

   double kinvis = 0.0005;
   double atm_pressure = 0.0;
    double gravity = 9.8;
    void init(const std::string &parameters_file);

} Parameters;

double get_next_parameter(std::ifstream &file);

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

void Acceleration_terms(const Vector &x, double t, Vector &a)
{
   double xi = x(0);
   double yi = x(1);
   double zi = x(2);

    a(0) = 0.0;
    a(1) = 0.0;
    a(2) = -Parameters.gravity;
}

int main(int argc, char *argv[])
{   
    //init Parameters
    std::string parameters_file = "settings/parameters.txt";
    Parameters.init(parameters_file);


    //Parameters
    double t = 0.0;
    int vis_print = 0;
    bool last_step = false;

    //Init MPI
    MPI_Session mpi(argc, argv);

    //Load Mesh
    Mesh mesh = Mesh("results/mesh.msh");
    mesh.EnsureNodes();
    int dim = mesh.Dimension();

    //Refine Serial Mesh
    for (int i = 0; i < Parameters.serial_refinements; ++i)
        mesh.UniformRefinement();

    //Make Parallel Mesh
    ParMesh pmesh = ParMesh(MPI_COMM_WORLD, mesh);
    
    //Refine Parallel Mesh
    for (int ii = 0; ii < Parameters.parallel_refinements; ii++)
        pmesh.UniformRefinement();

    //Create H1 Element Collection
    H1_FECollection vfec = H1_FECollection(Parameters.order, pmesh.Dimension());
    ParFiniteElementSpace vfes = ParFiniteElementSpace(&pmesh, &vfec, pmesh.Dimension());

    //Define Navier Solver
    NavierSolver flowsolver(&pmesh, Parameters.order, Parameters.kinvis);
    flowsolver.EnablePA(true);
    flowsolver.EnableNI(true);

    //Define Velocity Initial Conditions
    ParGridFunction *u_ic = flowsolver.GetCurrentVelocity();
    VectorFunctionCoefficient u_excoeff(pmesh.Dimension(), Initial_Velocity);
    u_ic->ProjectCoefficient(u_excoeff);

    //Define Pressure Initial Conditions
    ParGridFunction *u_ip = flowsolver.GetCurrentPressure();
    ConstantCoefficient  p_excoeff(Parameters.atm_pressure);
    u_ip->ProjectCoefficient(p_excoeff);

    //Define Velocity Boundary Conditions
    Array<int> attr(pmesh.bdr_attributes.Max());
    //inflow       sides           outflow
    attr[0] = 1;   attr[1] = 0;    attr[2] = 0;
    flowsolver.AddVelDirichletBC(Vel_Boundary_Condition, attr);

    //Define  Pressure Boundary Conditions
    Array<int> attr2(pmesh.bdr_attributes.Max());
    //inflow       sides           outflow
    attr2[0] = 0;  attr2[1] = 1;   attr2[2] = 1;
    flowsolver.AddPresDirichletBC(Press_Boundary_Condition, attr2);

    //Set Up Solver and Define Solution Pointers 
    flowsolver.Setup(Parameters.dt);
    ParGridFunction *u_gf = flowsolver.GetCurrentVelocity(); //Velocity solution
    ParGridFunction *p_gf = flowsolver.GetCurrentPressure(); //Presure solution

    //Calculate Vorticity
    ParGridFunction w_gf = ParGridFunction(&vfes);
    CurlGridFunctionCoefficient w(u_gf);
    w_gf.ProjectCoefficient(w);

    //Add Gravity Acceleration Term
    Array<int> domain_attr(pmesh.bdr_attributes.Max()); domain_attr=1;
    flowsolver.AddAccelTerm(Acceleration_terms,domain_attr);    

    //Paraview visualization
    ParaViewDataCollection paraview_out = ParaViewDataCollection("results/border_pressure", &pmesh);
    paraview_out.SetLevelsOfDetail(Parameters.order);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.RegisterField("Velocity", u_gf);
    paraview_out.RegisterField("Vorticity", &w_gf);
    paraview_out.RegisterField("Pressure", p_gf);
    paraview_out.SetCycle(vis_print);
    paraview_out.SetTime(t);
    paraview_out.Save();

    //Perform time integration
    for (int step = 0; !last_step; ++step)
    {
        //Check if last step
        if (t + Parameters.dt >= Parameters.t_final - Parameters.dt / 2)
            last_step = true;

        //Time step
        flowsolver.Step(t, Parameters.dt, step);

        //Update time in the Velocity Boundary condition
        u_excoeff.SetTime(t);

        //Compute CFL condition to check numerical stability
        //double cfl = flowsolver.ComputeCFL(*u_gf, Parameters.dt);
        //std::cout << step << "\t" << t << "\t" << Parameters.dt << "\t" << cfl << "\n";

        //Print Data for Visualization
        if (step%Parameters.vis_freq==0)
        {   
            vis_print++;
            CurlGridFunctionCoefficient w(u_gf);
            w_gf.ProjectCoefficient(w);
            paraview_out.SetCycle(vis_print);
            paraview_out.SetTime(t);
            paraview_out.Save();
        }
    }

    flowsolver.PrintTimingData();

    //with MPI sesion theres no neet to Finalize MPI

    //Free memory Dont have any pointer yet

    return 0;
}


double get_next_parameter(std::ifstream &file){
    std::string name_parameter;
    getline(file,name_parameter,' ');
    std::string string_value;
    getline(file,string_value);
    return stod(string_value);
}

void Navier_Parameters::init(const std::string &parameters_file){
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
