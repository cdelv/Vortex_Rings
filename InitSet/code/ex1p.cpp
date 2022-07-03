#include "mfem.hpp"
#include <fstream>
#include <string>
#include <iostream>

//using namespace std;
using namespace mfem;


struct Config
{
    //Numerical Method Parameters
    int n = 4;
    double Lx = 2, Ly = 3, Lz = 5;
    int serial_refinements = 1;
    int parallel_refinements = 0;
    int order = 2;
    int pid = -1;
    bool restart = false;
} Parameters;


void SaveState(ParMesh *pmesh,ParGridFunction *Velocity);
ParGridFunction InitState(ParMesh *pmesh,ParFiniteElementSpace &vfes);
ParMesh InitMesh();

double radius = 2;
double vel = 1;

void Initial_Velocity(const Vector &x, Vector &u){
  u(0) = x(0)*x(0);
  u(1) = 2*x(1);
  u(2) = 3*x(2);
}

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE.
   Mpi::Init();
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Parameters.pid = myid;
   Hypre::Init();

   
   // 2. Parse command-line options.
   int order = 2;
   bool static_cond = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = true;
   bool algebraic_ceed = false;

   // 3. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (Parameters.pid == 0) { device.Print(); }
   ParMesh pmesh_var = InitMesh();
   ParMesh *pmesh = &pmesh_var;

   return 0;
   H1_FECollection vfec = H1_FECollection(Parameters.order, pmesh->Dimension());
   ParFiniteElementSpace vfes = ParFiniteElementSpace(pmesh, &vfec, pmesh->Dimension());
   if (Parameters.pid == 0)
   {
     std::cout << "Number of finite element unknowns: "  << std::endl;
   }
   // 10. Define the solution vector x as a parallel finite element grid
   //     function corresponding to fespace. Initialize x with initial guess of
   //     zero, which satisfies the boundary conditions.
   ParGridFunction u_var = InitState(pmesh,vfes);
   ParGridFunction *u = &u_var;
   
   SaveState(pmesh,u);
   return 0;
}

//Print the final results
void SaveState(ParMesh *pmesh,ParGridFunction *Velocity){
    //Save final state
    std::ofstream out;
    std::ostringstream oss;
    oss << std::setw(10) << std::setfill('0') << Parameters.pid;
    std::string n_mesh = "Initial/pmesh_"+oss.str()+".msh";
    std::string n_velocity = "Initial/velocity"+oss.str()+".gf";
    
    out.open(n_mesh.c_str(),std::ios::out);
    pmesh->ParPrint(out);
    out.close();
    
    out.open(n_velocity.c_str(),std::ios::out);
    Velocity->Save(out);
    out.close();
}

ParGridFunction InitState(ParMesh *pmesh,ParFiniteElementSpace &vfes){
  ParGridFunction u(&vfes);
  if (Parameters.restart){
    VectorFunctionCoefficient u_init(pmesh->Dimension(), Initial_Velocity);
    u.ProjectCoefficient(u_init);
  }
  else{  
    std::ifstream in;
    std::ostringstream oss;
    oss << std::setw(10) << std::setfill('0') << Parameters.pid;
    std::string n_velocity = "Initial/velocity"+oss.str()+".gf";
    in.open(n_velocity.c_str(),std::ios::in);
    u = ParGridFunction(pmesh, in);
    in.close();
  }
  return u;
}

ParMesh InitMesh(){
  ParMesh pmesh;
  if(Parameters.restart){
    Mesh mesh = Mesh::MakeCartesian3D(2*Parameters.n, Parameters.n, Parameters.n, Element::QUADRILATERAL, Parameters.Lx, Parameters.Ly, Parameters.Lz);
    mesh.EnsureNodes();
    int dim = mesh.Dimension();

    //Refine Serial Mesh
    for (int i = 0; i < Parameters.serial_refinements; ++i)
      mesh.UniformRefinement();
    //Make Parallel Mesh
    //ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();
    //Refine mesh (parallel)
    for (int ii = 0; ii < Parameters.parallel_refinements; ii++)
      pmesh.UniformRefinement();
  } else {
    //Read the input mesh
    std::ifstream in;
    std::ostringstream oss;
    oss << std::setw(10) << std::setfill('0') << Parameters.pid;
    std::string n_mesh = "Initial/pmesh_"+oss.str()+".msh";
    in.open(n_mesh.c_str(),std::ios::in);
    pmesh.Load(in,1,0); 
    in.close();
  }
  return pmesh;
}
