#include "navier_solver.hpp"
#include <iostream>
using namespace mfem;
using namespace navier;

struct Navier_Parameters
{
    //Parameters
    int serial_refinements = 1;
    int parallel_refinements = 1;
    int order = 2;
    int vis_freq = 500;
    double dt = 0.0001;
    double t_final = 10.00;

    double R_ring = 0.2;
    double ring_center_y = 0.5;
    double ring_center_z = 0.5;
    double V0_ring = 0.0;
    double I_vel_exponent = 0.0;
    double I_active_dts = 100000.0;

   double kinvis = 0.0005;
   double atm_pressure = 0.0;
   double gravity = 9.8;

   double eta = 1e2;

} Parameters;

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

void gravity(const Vector &x, double t, Vector &a)
{
   double xi = x(0);
   double yi = x(1);
   double zi = x(2);

    a(0) = 0.0;
    a(1) = 0.0;
    a(2) = -Parameters.gravity;
}

class BrinkPenalAccel:public mfem::VectorCoefficient{
public:
    BrinkPenalAccel(int dim):mfem::VectorCoefficient(dim){}
    virtual ~BrinkPenalAccel(){}
    void SetVel(mfem::GridFunction* gfvel){vel=gfvel;}
    void SetTime(double tt){t=tt;}
    virtual void Eval(mfem::Vector &a, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
    {
    	//Body velocity  Fluid velocity   Fluid coordinates
        Vector U0;       Vector U;        Vector X;
        a.SetSize(GetVDim()); a*=0;
        //Set Body velocity
        U0.SetSize(GetVDim()); U0(0)=0; U0(1)=vy; U0(2)=vz;
        U.SetSize(GetVDim()); U*=0;

        if(t>10*Parameters.dt){
        	x=1.5+vx*(t-10*Parameters.dt);
        	U0(0)=vx;
        }

        //Get the physical coordinates of integration point
        T.Transform(ip,X); 

        if((std::pow(X(0)-x,2)+std::pow(X(1)-y,2)+std::pow(X(2)-z,2)<std::pow(r,2))){
            vel->GetVectorValue(T,ip,U);
            a=U;
            a-=U0;
            a*=-Parameters.eta;
        }
    }

private:
    mfem::GridFunction* vel = nullptr;
    double t = 0;
    double x  = 1.5;
    double y  = 0.5;
    double z  = 0.5;
    double r = 0.1;
    double vx  = -0.1;
    double vy  = 0.0;
    double vz  = 0.0;

};

int main(int argc, char *argv[])
{   
    //Parameters
    double t = 0.0;
    int vis_print = 0;
    bool last_step = false;

    //Init MPI
    MPI_Session mpi(argc, argv);

    //Load Mesh
    Mesh mesh = Mesh("mesh.msh");
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

    //Define Solution Pointers 
    ParGridFunction *u_gf = flowsolver.GetCurrentVelocity(); //Velocity solution
    ParGridFunction *p_gf = flowsolver.GetCurrentPressure(); //Presure solution

    //Calculate Vorticity
    ParGridFunction w_gf = ParGridFunction(&vfes);
    CurlGridFunctionCoefficient w(u_gf);
    w_gf.ProjectCoefficient(w);

    //Add Gravity Acceleration Term
    Array<int> domain_attr(pmesh.bdr_attributes.Max()); domain_attr=1;
    //flowsolver.AddAccelTerm(gravity,domain_attr);  

    //Create Solid object
    BrinkPenalAccel *piston= new BrinkPenalAccel(pmesh.Dimension());
    piston->SetVel(flowsolver.GetCurrentVelocity());
    flowsolver.AddAccelTerm(piston,domain_attr); 

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

    //Set Up solver
    flowsolver.Setup(Parameters.dt);

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

        piston->SetVel(flowsolver.GetCurrentVelocity());
        piston->SetTime(t);

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
    delete piston;

    return 0;
}
