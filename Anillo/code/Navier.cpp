#include "navier_solver.hpp"
#include <fstream>
using namespace mfem;
using namespace navier;

struct Navier_Parameters
{
    //Parameters
    int serial_refinements = 1;
    int parallel_refinements = 1;
    int order = 2;
    int vis_freq = 50;
    double dt = 0.001;
    double t_final = 20.00;

    double Ranillo = 0.3;
    double V0anillo = 7.0;

   double kinvis = 0.0005;
   double reference_pressure = 0.0;
   double Re = 1.0 / kinvis;
} Parameters;

void Initial_Velocity(const Vector &x, double t, Vector &u)
{
   double xi = x(0);
   double yi = x(1);
   double zi = x(2);

    u(0)=0;

  if(pow(yi- 0.5,2) + pow(zi - 0.5,2) < pow(Parameters.Ranillo,2)){
        u(0) = Parameters.V0anillo*exp(-40*xi);
  }

   u(1) = 0.0;
   u(2) = 0.0;
}

void Boundary_Condition1(const Vector &x, double t, Vector &u)
{
    double xi = x(0);
    double yi = x(1);
    double zi = x(2);

    if(pow(yi- 0.5,2) + pow(zi - 0.5,2) < pow(Parameters.Ranillo,2) && t<Parameters.dt*2){
        u(0) = Parameters.V0anillo;
    }
    else 
        u(0)=0;

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
    Mesh mesh = Mesh("mesh.msh");
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

    H1_FECollection vfec = H1_FECollection(Parameters.order, pmesh.Dimension());
    ParFiniteElementSpace vfes = ParFiniteElementSpace(&pmesh, &vfec, pmesh.Dimension());

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
    attr = 1; attr[0] = 1; attr[2] = 0;
    flowsolver.AddVelDirichletBC(Boundary_Condition1, attr);

    //Array<int> attr2(pmesh.bdr_attributes.Max());
    //attr2 = 1; attr2[0] = 0; attr2[1] = 0;
    //flowsolver.AddVelDirichletBC(Boundary_Condition2, attr2);

    flowsolver.Setup(Parameters.dt);

    ParGridFunction *u_gf = flowsolver.GetCurrentVelocity(); //Velocity solution
    ParGridFunction *p_gf = flowsolver.GetCurrentPressure(); //Presure solution

    //Vorticity
    ParGridFunction w_gf = ParGridFunction(&vfes);
    CurlGridFunctionCoefficient w(u_gf);
    w_gf.ProjectCoefficient(w);

    //Paraview visualization
    ParaViewDataCollection paraview_out = ParaViewDataCollection("results/graph", &pmesh);
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

        if (t + Parameters.dt >= Parameters.t_final - Parameters.dt / 2)
            last_step = true;

        flowsolver.Step(t, Parameters.dt, step);
        u_excoeff.SetTime(t);

        //double cfl = flowsolver.ComputeCFL(*u_gf, Parameters.dt);

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
