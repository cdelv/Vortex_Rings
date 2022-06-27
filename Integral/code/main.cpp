#include "navier_solver.hpp"
#include <fstream>
#include <string>
#include "boost/math/quadrature/gauss_kronrod.hpp"

using namespace mfem;
using namespace navier;
using namespace boost::math::quadrature;

//Configuration Functions
struct Config
{
    //Numerical Method Parameters
    int n = 10;
    int serial_refinements = 0;
    int parallel_refinements = 1;
    int order = 2;

    //Time Parameters
    int vis_freq = 100;
    double dt = 0.0001;
    double t_final = 1;

    //Box Parameters
    double Lx = 2;
    double Ly = 1;
    double Lz = 1;

    //Ring Parameters
    double R = 0.3;    //Radius
    double a = 0.1;   //Thickness
    double Rx = 1.;    //Position x
    double Ry = 0.5;   //Position y
    double Rz = 0.5;   //Position z
    double W = 100.;    //Mean Vorticity

    //Integral Parameters
    double Int_eps = 1E-9; 

    //Physical Parameters
    double kinvis = 1.48E-5;
    double atm_pressure = 0.;

    //Dimension Scale
    double CL = 1.;
    double CT = 1.;
    void Adimentionalize();

} Parameters;

//Assembly functions
void Initial_Vorticity(const Vector &x, double t, Vector &u);
void Initial_Velocity(const Vector &r, double t, Vector &u);
double integral(const Vector &r, double t, int coord);
void Vel_Boundary_Condition(const Vector &x, double t, Vector &u);
double Press_Boundary_Condition(const Vector &x, double t);

//Main function
int main(int argc, char *argv[])
{   
    //Init MPI
    MPI_Session mpi(argc, argv);

    //Aux Variables
    double t = 0.0;
    int vis_print = 0;
    bool last_step = false;
    //Parameters.Adimentionalize();
    NavierSolver *flowsolver = nullptr;

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

    //Calculate Vorticity
    ParGridFunction w = ParGridFunction(&vfes);
    VectorFunctionCoefficient w_bdr(pmesh->Dimension(), Initial_Vorticity);
    w.ProjectCoefficient(w_bdr);

    //Define Velocity Initial Conditions
    ParGridFunction *u = flowsolver->GetCurrentVelocity();
    VectorFunctionCoefficient u_init(pmesh->Dimension(), Initial_Velocity);
    u->ProjectCoefficient(u_init);

    //Define Pressure Initial Conditions
    ParGridFunction *p = flowsolver->GetCurrentPressure();
    ConstantCoefficient  p_init(Parameters.atm_pressure);
    p->ProjectCoefficient(p_init);

    //Create Domain Attr Array
    Array<int> domain_attr(pmesh->bdr_attributes.Max());
    domain_attr=1;

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

    //Set Up Solver (Must be After Adding Accel Terms)
    flowsolver->Setup(Parameters.dt);

    //Paraview Visualization
    ParaViewDataCollection paraview_out = ParaViewDataCollection("results/graph", pmesh);
    paraview_out.SetLevelsOfDetail(Parameters.order);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.RegisterField("Vorticity", &w);
    paraview_out.RegisterField("Pressure", p);
    paraview_out.RegisterField("Velocity", u);
    paraview_out.SetCycle(vis_print);
    paraview_out.SetTime(t);
    paraview_out.Save();

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

        if(mpi.Root()){
            std::cout.flush();
            std::cout << step << "\t" << Parameters.CT*t << "\t" << Parameters.dt << "\t" << vis_print << "\r";
        }

        //Print Data for Visualization
        if (step%Parameters.vis_freq==0 || last_step)
        {   
            vis_print++;
            CurlGridFunctionCoefficient u_curl(u);
            w.ProjectCoefficient(u_curl);
            paraview_out.SetCycle(vis_print);
            paraview_out.SetTime(Parameters.CT*t);
            paraview_out.Save();
        }
    }

    flowsolver->PrintTimingData();

    //Free memory
    delete pmesh;
    delete flowsolver;

    return 0;
}

void Config::Adimentionalize()
{
    CL = R;     //Critical Lenght
    CT = std::pow(2*R/a, 2)/(W*(std::log(8*R/4)-0.25)); //Critical Time

    dt /= CT;
    t_final /= CT;

    Lx /= CL;
    Ly /= CL;
    Lz /= CL;

    R  /= CL; 
    a  /= CL;
    Rx /= CL;
    Ry /= CL;
    Rz /= CL;
    W  *= CT;

    kinvis *= CT*pow(CL, -2);     
    atm_pressure *= pow(CT/CL, 2);
}

//Assembly functions
void LinealVortex(const double t, Vector &u){
    u(0) = Parameters.Rx;
    u(1) = Parameters.Ry + Parameters.R*std::cos(t);
    u(2) = Parameters.Rz + Parameters.R*std::sin(t);
}
 
void TangentVortex(const double t, Vector &u){
    u(0) = 0.;
    u(1) = -std::sin(t);
    u(2) = std::cos(t);
}
 
void Initial_Vorticity(const Vector &x, double t, Vector &u)
{ 
    double theta = std::atan2(x(2)-Parameters.Rz, x(1)-Parameters.Ry);
    LinealVortex(theta, u);   
    double r = u.DistanceTo(x);                        
 
    if (r <= Parameters.a){
        TangentVortex(theta, u);
        u *= (1-r/Parameters.a)*Parameters.W;   
    } else
        u = 0.;
}

void Initial_Velocity(const Vector &r, double t, Vector &u)
{
    u(0)=integral(r, t, 0);
    u(1)=integral(r, t, 1);
    u(2)=integral(r, t, 2);
}

double integral(const Vector &r, double t, int coord){
    double x1 = 0;
    double x2 = Parameters.Lx;
    double y1 = 0;
    double y2 = Parameters.Ly;
    double z1 = 0;
    double z2 = Parameters.Lz;

    const int points = 7;
    const int depth = 3;

    Vector cross; cross.SetSize(3);
    Vector rr; rr.SetSize(3);
    Vector W; W.SetSize(3);

    auto f2 = [&](double x, double y, double z){

        //W(r')
        rr(0)=x;
        rr(1)=y;
        rr(2)=z;
        Initial_Vorticity(rr,t,W);

        //r-r'
        rr(0)=r(0)-x;
        rr(1)=r(1)-y;
        rr(2)=r(2)-z;

        //w x (r-r')
        cross(0)= W(1)*rr(2)-W(2)*rr(1);
        cross(1)= W(2)*rr(0)-W(0)*rr(2);
        cross(2)= W(0)*rr(1)-W(1)*rr(0);

        double norm2 = rr.Norml2();

        if (norm2<1e-9)
            return 0.0;
        else
            return 1.0*cross(coord)/std::pow(norm2,1.5);
    };

    auto f1 = [&](double x, double y) { 

        auto g = [&](double z) { return f2(x,y,z); };

        return gauss_kronrod<double, points>::integrate(g, z1, z2, depth);
    };

    auto f = [&](double x) { 

        auto g = [&](double y) { return f1(x, y); };

        return gauss_kronrod<double, points>::integrate(g, y1, y2, depth);
    };

    double Q = gauss_kronrod<double, points>::integrate(f, x1, x2, depth);
    return 0.25*Q*M_1_PI;
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
