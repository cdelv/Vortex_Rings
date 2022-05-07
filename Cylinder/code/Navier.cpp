#include "navier_solver.hpp"
#include <fstream>
using namespace mfem;
using namespace navier;

struct Navier_Parameters
{
    //Parameters
    int serial_refinements = 1;
    int parallel_refinements = 1;
    int order = 4;
    int vis_freq = 500;
    double dt = 0.0001;
    double t_final = 8.0;

   double kinvis = 0.001;
   double reference_pressure = 0.0;
   double Re = 1.0 / kinvis;
   double lam = 0.5*Re-sqrt(0.25*Re*Re+4.0*M_PI*M_PI);
} Parameters;

void Initial_Velocity(const Vector &x, double t, Vector &u)
{
   double xi = x(0);
   double yi = x(1);
   double zi = x(2);

   double U = 2.25;

   if (xi <= 1e-8)
   {
      u(0) = 16.0 * U * yi * zi * sin(M_PI * t / 8.0) * (0.41 - yi)
             * (0.41 - zi) / pow(0.41, 4.0);
   }
   else
   {
      u(0) = 0.0;
   }
   u(1) = 0.0;
   u(2) = 0.0;
}

double Initial_Pressure(const Vector &x, double t)
{
    double xi = x(0);
    double yi = x(1);
    double zi = x(2);

    return 0.00001;
   //return 0.5*(1.0-exp(2.0*Parameters.lam*xi))+Parameters.reference_pressure;
}

void Boundary_Condition(const Vector &x, double t, Vector &u)
{
    double xi = x(0);
    double yi = x(1);
    double zi = x(2);


    u(0) = (t*t+0.01)*0.005/(t+1);
    u(1) = 0.0;
    u(2) = 0.0;
}

int main(int argc, char *argv[])
{   
    //Parameters
    double t = 0.0;
    int vis_print = 0;
    bool last_step = false;

    //Init MPI
    MPI_Session mpi(argc, argv);

    //Load mesh
    //Mesh mfem::Mesh::MakeCartesian3D(int nx,int ny,int nz,Element::Type type,double sx=1.0,double sy=1.0,double sz=1.0, bool sfc_ordering=true)   
    //Mesh mesh = Mesh::MakeCartesian3D(2,4,2,Element::QUADRILATERAL,1.5,2.0,0.1);
    Mesh mesh = Mesh("box-cylinder.mesh");
    mesh.EnsureNodes();
    int dim = mesh.Dimension();

    //Refine serial mesh
    for (int i = 0; i < Parameters.serial_refinements; ++i)
        mesh.UniformRefinement();

    //Make parallel mesh
    ParMesh pmesh = ParMesh(MPI_COMM_WORLD, mesh);
    
    //Refine parallel mesh
    for (int ii = 0; ii < Parameters.parallel_refinements; ii++)
        pmesh.UniformRefinement();

    //Define NavierSolver
    NavierSolver flowsolver(&pmesh, Parameters.order, Parameters.kinvis);
    flowsolver.EnablePA(true);
    flowsolver.EnableNI(true);

    //Define Initial conditions
    ParGridFunction *u_ic = flowsolver.GetCurrentVelocity();
    VectorFunctionCoefficient u_excoeff(pmesh.Dimension(), Initial_Velocity);
    u_ic->ProjectCoefficient(u_excoeff);

    //Pressure
    //FunctionCoefficient p_excoeff(Initial_Pressure);

    //Define Dirichlet boundary conditions
    Array<int> attr(pmesh.bdr_attributes.Max());
    attr[0] = 1;
    attr[2] = 1;
    flowsolver.AddVelDirichletBC(Initial_Velocity, attr);

     flowsolver.Setup(Parameters.dt);

    ParGridFunction *u_gf = flowsolver.GetCurrentVelocity(); //Velocity solution
    ParGridFunction *p_gf = flowsolver.GetCurrentPressure(); //Presure solution

    //Paraview visualization
    ParaViewDataCollection paraview_out = ParaViewDataCollection("results/graph", &pmesh);
    paraview_out.SetLevelsOfDetail(Parameters.order);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.RegisterField("Velocity", u_gf);
    paraview_out.RegisterField("Pressure", p_gf);
    paraview_out.SetCycle(vis_print);
    paraview_out.SetTime(t);
    paraview_out.Save();

    //Perform time integration
    for (int step = 0; !last_step; ++step)
    {

        if (t + Parameters.dt >= Parameters.t_final - Parameters.dt / 2)
            last_step = true;

        flowsolver.Step(t, Parameters.dt, step);

        //double cfl = flowsolver.ComputeCFL(*u_gf, Parameters.dt);

        if (step%Parameters.vis_freq==0)
        {   
            vis_print++;
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