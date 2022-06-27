#include "navier_solver.hpp"
#include <fstream>
#include <string>
#include "boost/math/quadrature/gauss_kronrod.hpp"

using namespace mfem;
using namespace navier;
using namespace boost::math::quadrature;

//Number of Integration Points
const int points = 7;

//Configuration Functions
struct Config
{
    //Numerical Method Parameters
    int n = 6;
    int serial_refinements = 1;
    int parallel_refinements = 1;
    int order = 2;

    //Time Parameters
    int vis_freq = 1000;
    double dt = 0.0001;
    double t_final = 3;

    //Box Parameters
    double Lx = 4.0;
    double Ly = 1.5;
    double Lz = 1.5;

    //Ring Parameters
    double R = 0.3;          //Radius
    double a = 0.1;         //Thickness
    double Rx = Lx*0.3;    //Position x
    double Ry = Ly*0.5;   //Position y
    double Rz = Lz*0.5;  //Position z
    double W = 100.;    //Mean Vorticity

    //Integral Parameters
    double Int_eps = 1E-15; 
    double Int_cutoff = 1.0;
    int depth = 3;

    //Physical Parameters
    double kinvis = 1.48E-5;
    double atm_pressure = 0.;

    //Dimension Scale
    double CL = R;                                              //Critical Lenght
    double CT = std::pow(2*R/a, 2)/(W*(std::log(8*R/a)-0.25)); //Critical Time
    void Adimentionalize();
} Parameters;

//Assembly Functions
void Initial_Vorticity(const Vector &x, double t, Vector &u);
void Initial_Velocity(const Vector &r, double t, Vector &u);
double Integral(const Vector &r, double t, int coord);
void Vel_Boundary_Condition(const Vector &x, double t, Vector &u);
double Press_Boundary_Condition(const Vector &x, double t);

//Main Function
int main(int argc, char *argv[])
{   
    //Init MPI
    MPI_Session mpi(argc, argv);

    //Aux Variables
    double t = 0.0;
    int vis_print = 0;
    bool last_step = false;
    NavierSolver *flowsolver = nullptr;

    //Adimentionalize
    Parameters.Adimentionalize();

    ParMesh *pmesh = new ParMesh();
    {
    //Load Mesh (In Different Scope to Delete it After Parallel Mesh is Created)
    Mesh mesh = Mesh::MakeCartesian3D(2*Parameters.n, Parameters.n, Parameters.n, Element::QUADRILATERAL, Parameters.Lx, Parameters.Ly, Parameters.Lz);
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
    ParFiniteElementSpace pfes = ParFiniteElementSpace(pmesh, &pfec);                    //Scalar Space

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
    attr[4] = 0;  attr[2] = 0;   attr[0] = 0;   attr[5] = 0;  attr[1] = 0;   attr[3] = 0;
    flowsolver->AddVelDirichletBC(Vel_Boundary_Condition, attr);

    //Define  Pressure Boundary Conditions
    Array<int> attr2(pmesh->bdr_attributes.Max());
    //Inflow       Outflow         Side down       Side up        Side left       Side right
    attr2[4] = 1;  attr2[2] = 1;   attr2[0] = 1;   attr2[5] = 1;  attr2[1] = 1;   attr2[3] = 1;
    flowsolver->AddPresDirichletBC(Press_Boundary_Condition, attr2);

    //Set Up Solver (After Adding Accel Terms)
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

    //Free Memory
    delete pmesh;
    delete flowsolver;

    return 0;
}

void Config::Adimentionalize()
{
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

    Int_cutoff /= CL; 

    kinvis *= CT*pow(CL, -2);     
    atm_pressure *= pow(CT/CL, 2);
}

//Assembly Functions
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
    u(0)=Integral(r, t, 0);
    u(1)=Integral(r, t, 1);
    u(2)=Integral(r, t, 2);
}

double Integral(const Vector &r, double t, int coord){

    Vector W; W.SetSize(3);
    double theta = std::atan2(r(2)-Parameters.Rz, r(1)-Parameters.Ry);
    LinealVortex(theta, W);   
    double d = W.DistanceTo(r); 
    if (d>Parameters.Int_cutoff)
        return 0.0;

    double x1 = Parameters.Rx-1.5*Parameters.a;
    double x2 = Parameters.Rx+1.5*Parameters.a;
    double y1 = Parameters.Ry-Parameters.R-1.5*Parameters.a;
    double y2 = Parameters.Ry+Parameters.R+1.5*Parameters.a;
    double z1 = Parameters.Rz-Parameters.R-1.5*Parameters.a;
    double z2 = Parameters.Rz+Parameters.R+1.5*Parameters.a;

    Vector cross; cross.SetSize(3);
    Vector rr; rr.SetSize(3);

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

        //W(r') x (r-r')
        cross(0)= W(1)*rr(2)-W(2)*rr(1);
        cross(1)= W(2)*rr(0)-W(0)*rr(2);
        cross(2)= W(0)*rr(1)-W(1)*rr(0);

        double norm = std::sqrt(rr(0)*rr(0)+rr(1)*rr(1)+rr(2)*rr(2));

        if (norm<Parameters.Int_eps)
            return 0.0;
        else
            return cross(coord)/std::pow(norm,3);
    };

    auto f1 = [&](double x, double y) { 

        auto g = [&](double z) { return f2(x,y,z); };

        return gauss_kronrod<double, points>::integrate(g, z1, z2, Parameters.depth);
    };

    auto f = [&](double x) { 

        auto g = [&](double y) { return f1(x, y); };

        return gauss_kronrod<double, points>::integrate(g, y1, y2, Parameters.depth);
    };

    return 0.25*M_1_PI*gauss_kronrod<double, points>::integrate(f, x1, x2, Parameters.depth);
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
