#include "InitSet.h"

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

void InitState(ParMesh *pmesh,ParGridFunction *u){
  std::ifstream in;
  std::ostringstream oss;
  oss << std::setw(10) << std::setfill('0') << Parameters.pid;
  std::string n_velocity = "Initial/velocity"+oss.str()+".gf";
  in.open(n_velocity.c_str(),std::ios::in);
  u = new ParGridFunction(pmesh, in);
  in.close();
}

void InitMesh(){
  Mesh mesh = Mesh::MakeCartesian3D(2*Parameters.n, Parameters.n, Parameters.n, Element::QUADRILATERAL, Parameters.Lx, Parameters.Ly, Parameters.Lz);
  mesh.EnsureNodes();
  int dim = mesh.Dimension();

  //Refine Serial Mesh
  for (int i = 0; i < Parameters.serial_refinements; ++i)
    mesh.UniformRefinement();

  //Make Parallel Mesh
  pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
  mesh.Clear();
  if (!Parameters.restart){
    //Refine mesh (parallel)
    for (int ii = 0; ii < Parameters.parallel_refinements; ii++)
      pmesh->UniformRefinement();
  } else {
    //Read the input mesh
    std::ifstream in;
    std::ostringstream oss;
    oss << std::setw(10) << std::setfill('0') << Parameters.pid;
    std::string n_mesh = "Initial/pmesh_"+oss.str()+".msh";
    in.open(n_mesh.c_str(),std::ios::in);
    pmesh->Load(in,1,0); 
    in.close();
  }
}
