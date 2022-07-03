#include <fstream>
#include <string>
#include "../Navier/navier_solver.hpp"
#include <boost/math/tools/minima.hpp>
#include "boost/math/quadrature/gauss_kronrod.hpp"

using namespace mfem;
using namespace navier;
using namespace boost::math::quadrature;

//Number of Integration Points
const int points = 25;

//Configuration Functions
struct Config
{
    //Numerical Method Parameters
    int n = 4;
    int serial_refinements = 0;
    int parallel_refinements = 1;
    int order = 2;

    //Time Parameters
    int vis_freq = 1000;
    double dt = 0.0001;
    double t_final = 16;

    //Box Parameters
    double Lx = 4.0;
    double Ly = 1.5;
    double Lz = 1.5;

    //Ring Parameters
    double R = 0.3;           //Radius
    double epsilon = 0.1;    //Perturbation
    double a = 0.1;         //Thickness
    double Rx = Lx*0.3;    //Position x
    double Ry = Ly*0.5;   //Position y
    double Rz = Lz*0.5;  //Position z
    double W = 10.0;    //Mean Vorticity

    //Integral Parameters
    double Int_eps = 1E-14; 
    int depth = 3;

    //Physical Parameters
    double kinvis = 1.48E-4;
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
void Compute_Curl_Error(ParMesh *pmesh, ParGridFunction *u, ParGridFunction w, VectorFunctionCoefficient w_bdr, bool print);

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
    ParaViewDataCollection paraview_out = ParaViewDataCollection("results/ellipse", pmesh);
    paraview_out.SetLevelsOfDetail(Parameters.order);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.RegisterField("Vorticity", &w);
    paraview_out.RegisterField("Pressure", p);
    paraview_out.RegisterField("Velocity", u);
    paraview_out.SetCycle(vis_print);
    paraview_out.SetTime(t);
    paraview_out.Save();

    Compute_Curl_Error(pmesh, u, w, w_bdr,mpi.Root());

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
    epsilon /= CL; 
    a  /= CL;
    Rx /= CL;
    Ry /= CL;
    Rz /= CL;
    W  *= CT;

    kinvis *= CT*pow(CL, -2);     
    atm_pressure *= pow(CT/CL, 2);
}

//Assembly Functions
void LinealVortex(const double t, Vector &u){
    u(0) = Parameters.Rx;
    u(1) = Parameters.Ry + (Parameters.R + Parameters.epsilon)*std::cos(t);
    u(2) = Parameters.Rz + (Parameters.R - Parameters.epsilon)*std::sin(t);
}
 
void TangentVortex(const double t, Vector &u){
    u(0) = 0.;
    u(1) = -(Parameters.R + Parameters.epsilon)*std::sin(t);
    u(2) = (Parameters.R - Parameters.epsilon)*std::cos(t);
    double norm = u.Norml2();
    u /= norm;
}
 
void Initial_Vorticity(const Vector &x, double t, Vector &u)
{ 
    auto f = [&](double theta){
        LinealVortex(theta, u);   
        return u.DistanceTo(x);                        
    };

    double theta = std::atan2(x(2)-Parameters.Rz, x(1)-Parameters.Ry);
    double r = 0.;

    int bits = std::numeric_limits<double>::digits;
    std::pair<double, double> mini = boost::math::tools::brent_find_minima(f, theta-M_PI_4, theta+M_PI_4, bits);

    theta = mini.first;
    r = mini.second;
 
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

        double norm = rr.Norml2();
        return cross(coord)*std::pow(norm+Parameters.Int_eps,-3);
    };

    auto f1 = [&](double x, double y) { 

        auto g = [&](double z) { return f2(x,y,z); };

        //return gauss_kronrod<double, points>::integrate(g, z1, z2, Parameters.depth);
        return gauss<double, points>::integrate(g, z1, z2);
    };

    auto f = [&](double x) { 

        auto g = [&](double y) { return f1(x, y); };

        //return gauss_kronrod<double, points>::integrate(g, y1, y2, Parameters.depth);
        return gauss<double, points>::integrate(g, y1, y2);
    };

    //return 0.25*M_1_PI*gauss_kronrod<double, points-4>::integrate(f, x1, x2, Parameters.depth);
    return 0.25*M_1_PI*gauss<double, points-18>::integrate(f, x1, x2);
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

void Compute_Curl_Error(ParMesh *pmesh, ParGridFunction *u, ParGridFunction w, VectorFunctionCoefficient w_bdr, bool print)
{
    //Create Integration Rule
    const IntegrationRule *irs[Geometry::NumGeom];
    for (int i=0; i < Geometry::NumGeom; ++i)
        irs[i] = &(IntRules.Get(i, 2*Parameters.order-1));

    //Compute L2 Norm of Initial VorticityCoe
    double norm = ComputeGlobalLpNorm(2,w_bdr,*pmesh,irs);

    //Compute Velocity Curl
    CurlGridFunctionCoefficient u_curl(u);

    //Compute Error
    double Error = w.ComputeL2Error(u_curl,irs)/norm;

    if(print){
        std::cout <<"Initial Velocity Curl L2 Relative Error = "<< Error << std::endl;
        std::cout <<"L2 Vorticity Norm = "<< norm << std::endl;
    }
}