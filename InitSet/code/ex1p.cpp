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

  
  //We have to include this two variables
    int pid = -1;
    bool restart = true;

} Parameters;


void SaveState(ParMesh *pmesh,ParGridFunction *Velocity);
ParGridFunction InitState(ParMesh *pmesh,ParFiniteElementSpace &vfes);
ParMesh InitMesh();

void Initial_Velocity(const Vector &x, Vector &u){
  u(0) = x(0)*x(0);
  u(1) = 2*x(1);
  u(2) = 3*x(2);
}

int main(int argc, char *argv[])
{
   // Initialize MPI and HYPRE.
   Mpi::Init();
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Parameters.pid = myid;
   Hypre::Init();

   
   // Parse command-line options.
   int order = 2;
   bool static_cond = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = true;
   bool algebraic_ceed = false;
   Device device(device_config);
   if (Parameters.pid == 0) { device.Print(); }


   //Parallel MESH
   ParMesh pmesh_var = InitMesh();
   ParMesh *pmesh = &pmesh_var;
   H1_FECollection vfec = H1_FECollection(Parameters.order, pmesh->Dimension());
   ParFiniteElementSpace vfes = ParFiniteElementSpace(pmesh, &vfec, pmesh->Dimension());
   
   
   //Define velocity coeff in order to checkcomparison the error with the coeff readed from file (when paramaters.restart = false)
   VectorFunctionCoefficient u_check(pmesh->Dimension(), Initial_Velocity);
   //Create or read the initial velocity
   ParGridFunction u_var = InitState(pmesh,vfes);
   ParGridFunction *u = &u_var;

   //get error
   double error = u_var.ComputeL2Error(u_check);
   //print error 
   if(myid == 0){
     std::cout << "Restart Flag:\t" << Parameters.restart << "\n";
     std::cout << "Error\t" << error << "\n";
   }

   //save the state if restart = true 
   SaveState(pmesh,u);
   
   return 0;
}

//Print the final results
void SaveState(ParMesh *pmesh,ParGridFunction *Velocity){

  if(Parameters.restart){
    //Save final state
    std::ofstream out;
    std::ostringstream oss;
    oss << std::setw(5) << std::setfill('0') << Parameters.pid;
    std::string n_mesh = "Initial/pmesh_"+oss.str()+".msh";
    const char *n_velocity = "Initial/velocity";
    
    out.open(n_mesh.c_str(),std::ios::out);
    out.precision(16);
    pmesh->ParPrint(out);
    out.close();
    
    Velocity->Save(n_velocity);
  }
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
    oss << Parameters.pid;
    std::string n_velocity = "Initial/velocity.00000"+std::to_string(Parameters.pid);
    in.open(n_velocity.c_str(),std::ios::in);
    u = ParGridFunction(pmesh, in);
    in.close();
  }
  return u;
}

ParMesh InitMesh(){

  //if there are not files yet
  if(Parameters.restart){
    Mesh mesh = Mesh::MakeCartesian3D(2*Parameters.n, Parameters.n, Parameters.n, Element::QUADRILATERAL, Parameters.Lx, Parameters.Ly, Parameters.Lz);
    mesh.EnsureNodes();
    int dim = mesh.Dimension();

    //Refine Serial Mesh
    for (int i = 0; i < Parameters.serial_refinements; ++i)
      mesh.UniformRefinement();
    //Make Parallel Mesh
    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();
    //Refine mesh (parallel)
    for (int ii = 0; ii < Parameters.parallel_refinements; ii++)
      pmesh.UniformRefinement();
    return pmesh;
    } else {
    //Read the input mesh
    std::ostringstream oss;
    oss << std::setw(5) << std::setfill('0') << Parameters.pid;
    std::string n_mesh = "Initial/pmesh_"+oss.str()+".msh";
    std::ifstream mesh_ifs(n_mesh.c_str());
    ParMesh pmesh(MPI_COMM_WORLD, mesh_ifs);
    return pmesh;
    }
}
