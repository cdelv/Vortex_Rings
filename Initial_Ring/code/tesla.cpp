// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.
//
//            -----------------------------------------------------
//            Tesla Miniapp:  Simple Magnetostatics Simulation Code
//            -----------------------------------------------------
//
// This miniapp solves a simple 3D magnetostatic problem.
//
//                     Curl 1/mu Curl A = J + Curl mu0/mu M
//
// The permeability function is that of the vacuum with an optional diamagnetic
// or paramagnetic spherical shell. The optional current density takes the form
// of a user defined ring of current. The optional magnetization consists of a
// cylindrical bar of constant magnetization.
//
// The boundary conditions either apply a user selected uniform magnetic flux
// density or a surface current flowing between user defined surfaces.
//
// We discretize the vector potential with H(Curl) finite elements. The magnetic
// flux B is discretized with H(Div) finite elements.
//
// Compile with: make tesla
//
// Sample runs:
//
//   A cylindrical bar magnet in a metal sphere:
//      mpirun -np 4 tesla -bm '0 -0.5 0 0 0.5 0 0.2 1'
//
//   A spherical shell of paramagnetic material in a uniform B field:
//      mpirun -np 4 tesla -ubbc '0 0 1' -ms '0 0 0 0.2 0.4 10'
//
//   A ring of current in a metal sphere:
//      mpirun -np 4 tesla -cr '0 0 -0.2 0 0 0.2 0.2 0.4 1'
//
//   A Halbach array of permanent magnets:
//      mpirun -np 4 tesla -m ../../data/beam-hex.mesh -rs 2
//                         -ha '1 0.1 0.3 7 0.9 0.7 0 1 12'
//
//   An example demonstrating the use of surface currents:
//      mpirun -np 4 tesla -m square-angled-pipe.mesh
//                         -kbcs '3' -vbcs '1 2' -vbcv '-0.5 0.5'
//
//   An example combining the paramagnetic shell, permanent magnet,
//   and current ring:
//      mpirun -np 4 tesla -m ../../data/inline-hex.mesh
//                         -ms '0.5 0.5 0.5 0.4 0.45 20'
//                         -bm '0.5 0.5 0.3 0.5 0.5 0.7 0.1 1'
//                         -cr '0.5 0.5 0.45 0.5 0.5 0.55 0.2 0.3 1'
//
//   By default the sources and fields are all zero:
//      mpirun -np 4 tesla

#include "tesla_solver.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;
using namespace mfem::electromagnetics;

void current_ring(const Vector &, Vector &);

int main(int argc, char *argv[])
{
    MPI_Session mpi(argc, argv);

   // Parse command-line options.
   int order = 1;
   int serial_ref_levels = 1;
   int parallel_ref_levels = 2;
   bool visualization = true;
   bool visit = true;

   Array<int> kbcs;
   Array<int> vbcs;
   Vector vbcv;

   ParMesh *pmesh = new ParMesh();
   {
   //Load Mesh (Pointer to Delete it After Parallel Mesh is Created)
   int n = 10;
   double Lx = 5.;   
   double Ly = 5.;
   double Lz = 5.;
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::QUADRILATERAL, Lx, Ly, Lz);
   mesh.EnsureNodes();
   int dim = mesh.Dimension();

   //Refine Serial Mesh
   for (int i = 0; i < serial_ref_levels; ++i)
       mesh.UniformRefinement();

   //Make Parallel Mesh
   pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
   }
   
   //Refine Parallel Mesh
   for (int ii = 0; ii < parallel_ref_levels; ii++)
       pmesh->UniformRefinement();

   H1_FECollection vfec = H1_FECollection(order, pmesh->Dimension());
   ParFiniteElementSpace vfes = ParFiniteElementSpace(pmesh, &vfec, pmesh->Dimension());
   ParGridFunction w = ParGridFunction(&vfes);
   VectorFunctionCoefficient w_bdr(pmesh->Dimension(), current_ring);
   w.ProjectCoefficient(w_bdr);

   // Create a coefficient describing the magnetic permeability
   ConstantCoefficient muInvCoef(1.);

   // Create the Magnetostatic solver
   TeslaSolver Tesla(*pmesh, order, kbcs, vbcs, vbcv, muInvCoef, NULL, current_ring, NULL);

   // Display the current number of DoFs in each finite element space
   Tesla.PrintSizes();

   // Assemble all forms
   Tesla.Assemble();

   // Solve the system and compute any auxiliary fields
   Tesla.Solve();

   ParGridFunction psi(Tesla.GetVectorPotential());
   ParGridFunction v(Tesla.GetVectorField());

   ParaViewDataCollection paraview_out = ParaViewDataCollection("results/graph", pmesh);
   paraview_out.SetLevelsOfDetail(order);
   paraview_out.SetDataFormat(VTKFormat::BINARY);
   paraview_out.RegisterField("Vorticity", &w);
   paraview_out.RegisterField("Psi", &psi);
   paraview_out.RegisterField("Velocity", &v);
   paraview_out.SetCycle(0);
   paraview_out.SetTime(0);
   paraview_out.Save();

   // Determine the current size of the linear system
   int prob_size = Tesla.GetProblemSize();

   return 0;
}

// An annular ring of current density.  The ring has two axis end
// points, inner and outer radii, and a constant current in Amperes.
void current_ring(const Vector &x, Vector &j)
{

   // Current Density Function
   // Axis Start, Axis End, Inner Ring Radius,
   //                               Outer Ring Radius, and Total Current
   //                               of current ring (annulus)
   Vector cr_params_(9);
   //      mpirun -np 4 tesla -cr '0 0 -0.2 0 0 0.2 0.2 0.4 1'
   cr_params_[0] = 0; cr_params_[1] = 2.5; cr_params_[2] = 2.5; 
   cr_params_[3] = 5; cr_params_[4] = 2.5; cr_params_[5] = 2.5; 
   
   cr_params_[6] = 1; cr_params_[7] = 2; 

   cr_params_[8] = 1.; 

   j.SetSize(x.Size());
   j = 0.0;

   Vector  a(x.Size());  // Normalized Axis vector
   Vector xu(x.Size());  // x vector relative to the axis end-point
   Vector ju(x.Size());  // Unit vector in direction of current

   xu = x;

   for (int i=0; i<x.Size(); i++)
   {
      xu[i] -= cr_params_[i];
      a[i]   = cr_params_[x.Size()+i] - cr_params_[i];
   }

   double h = a.Norml2();

   if ( h == 0.0 )
   {
      return;
   }

   double ra = cr_params_[2*x.Size()+0];
   double rb = cr_params_[2*x.Size()+1];
   if ( ra > rb )
   {
      double rc = ra;
      ra = rb;
      rb = rc;
   }
   double xa = xu*a;

   if ( h > 0.0 )
   {
      xu.Add(-xa/(h*h),a);
   }

   double xp = xu.Norml2();

   if ( xa >= 0.0 && xa <= h*h && xp >= ra && xp <= rb )
   {
      ju(0) = a(1) * xu(2) - a(2) * xu(1);
      ju(1) = a(2) * xu(0) - a(0) * xu(2);
      ju(2) = a(0) * xu(1) - a(1) * xu(0);
      ju /= h;

      j.Add(cr_params_[2*x.Size()+2]/(h*(rb-ra)),ju);
   }
}
