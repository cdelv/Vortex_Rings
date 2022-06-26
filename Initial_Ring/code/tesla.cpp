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

// Permeability Function
Coefficient * SetupInvPermeabilityCoefficient();

static Vector pw_mu_(0);      // Piecewise permeability values
static Vector pw_mu_inv_(0);  // Piecewise inverse permeability values
static Vector ms_params_(0);  // Center, Inner and Outer Radii, and
//                               Permeability of magnetic shell
double magnetic_shell(const Vector &);
double magnetic_shell_inv(const Vector & x) { return 1.0/magnetic_shell(x); }

void current_ring(const Vector &, Vector &);

// Magnetization
static Vector bm_params_(0);  // Axis Start, Axis End, Bar Radius,
//                               and Magnetic Field Magnitude
void bar_magnet(const Vector &, Vector &);

static Vector ha_params_(0);  // Bounding box,
//                               axis index (0->'x', 1->'y', 2->'z'),
//                               rotation axis index
//                               and number of segments
void halbach_array(const Vector &, Vector &);

// A Field Boundary Condition for B = (Bx,By,Bz)
static Vector b_uniform_(0);
void a_bc_uniform(const Vector &, Vector&);

// Phi_M Boundary Condition for H = (0,0,1)
double phi_m_bc_uniform(const Vector &x);

// Prints the program's logo to the given output stream
void display_banner(ostream & os);

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();

   if ( Mpi::Root() ) { display_banner(cout); }

   // Parse command-line options.
   const char *mesh_file = "../../data/ball-nurbs.mesh";
   int order = 1;
   int maxit = 100;
   int serial_ref_levels = 0;
   int parallel_ref_levels = 0;
   bool visualization = true;
   bool visit = true;

   Array<int> kbcs;
   Array<int> vbcs;

   Vector vbcv;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&serial_ref_levels, "-rs", "--serial-ref-levels",
                  "Number of serial refinement levels.");
   args.AddOption(&parallel_ref_levels, "-rp", "--parallel-ref-levels",
                  "Number of parallel refinement levels.");
   args.AddOption(&b_uniform_, "-ubbc", "--uniform-b-bc",
                  "Specify if the three components of the constant magnetic flux density");
   args.AddOption(&pw_mu_, "-pwm", "--piecewise-mu",
                  "Piecewise values of Permeability");
   args.AddOption(&ms_params_, "-ms", "--magnetic-shell-params",
                  "Center, Inner Radius, Outer Radius, and Permeability of Magnetic Shell");
   //args.AddOption(&cr_params_, "-cr", "--current-ring-params",
   //               "Axis End Points, Inner Radius, Outer Radius and Total Current of Annulus");
   args.AddOption(&bm_params_, "-bm", "--bar-magnet-params",
                  "Axis End Points, Radius, and Magnetic Field of Cylindrical Magnet");
   args.AddOption(&ha_params_, "-ha", "--halbach-array-params",
                  "Bounding Box Corners and Number of Segments");
   args.AddOption(&kbcs, "-kbcs", "--surface-current-bc",
                  "Surfaces for the Surface Current (K) Boundary Condition");
   args.AddOption(&vbcs, "-vbcs", "--voltage-bc-surf",
                  "Voltage Boundary Condition Surfaces (to drive K)");
   args.AddOption(&vbcv, "-vbcv", "--voltage-bc-vals",
                  "Voltage Boundary Condition Values (to drive K)");
   args.AddOption(&maxit, "-maxit", "--max-amr-iterations",
                  "Max number of iterations in the main AMR loop.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&visit, "-visit", "--visit", "-no-visit", "--no-visit",
                  "Enable or disable VisIt visualization.");
   args.Parse();
   if (!args.Good())
   {
      if (Mpi::Root())
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (Mpi::Root())
   {
      args.PrintOptions(cout);
   }

   ParMesh *pmesh = new ParMesh();
   {
   //Load Mesh (Pointer to Delete it After Parallel Mesh is Created)
   int n = 6;
   double Lx = 2.;   
   double Ly = 1.;
   double Lz = 1.;
   Mesh mesh = Mesh::MakeCartesian3D(2*n, n, n, Element::QUADRILATERAL, Lx, Ly, Lz);
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

   // Create a coefficient describing the magnetic permeability
   ConstantCoefficient muInvCoef(1.);

   // Create the Magnetostatic solver
   TeslaSolver Tesla(*pmesh, order, kbcs, vbcs, vbcv, muInvCoef,
           NULL, current_ring, NULL);
                 /*    (b_uniform_.Size() > 0 ) ? a_bc_uniform  : NULL,
                     (cr_params_.Size() > 0 ) ? current_ring  : NULL,
                     (bm_params_.Size() > 0 ) ? bar_magnet    : 
                     (ha_params_.Size() > 0 ) ? halbach_array : NULL);*/


   // Display the current number of DoFs in each finite element space
   Tesla.PrintSizes();

   // Assemble all forms
   Tesla.Assemble();

   // Solve the system and compute any auxiliary fields
   Tesla.Solve();

   // Determine the current size of the linear system
   int prob_size = Tesla.GetProblemSize();

   return 0;
}

// Print the Volta ascii logo to the given ostream
void display_banner(ostream & os)
{
   os << "  ___________            __            " << endl
      << "  \\__    ___/___   _____|  | _____     " << endl
      << "    |    |_/ __ \\ /  ___/  | \\__  \\    " << endl
      << "    |    |\\  ___/ \\___ \\|  |__/ __ \\_  " << endl
      << "    |____| \\___  >____  >____(____  /  " << endl
      << "               \\/     \\/          \\/   " << endl << flush;
}

// The Permeability is a required coefficient which may be defined in
// various ways so we'll determine the appropriate coefficient type here.
Coefficient *
SetupInvPermeabilityCoefficient()
{
   Coefficient * coef = NULL;

   if ( ms_params_.Size() > 0 )
   {
      coef = new FunctionCoefficient(magnetic_shell_inv);
   }
   else if ( pw_mu_.Size() > 0 )
   {
      pw_mu_inv_.SetSize(pw_mu_.Size());
      for (int i = 0; i < pw_mu_.Size(); i++)
      {
         MFEM_ASSERT( pw_mu_[i] > 0.0, "permeability values must be positive" );
         pw_mu_inv_[i] = 1.0/pw_mu_[i];
      }
      coef = new PWConstCoefficient(pw_mu_inv_);
   }
   else
   {
      coef = new ConstantCoefficient(1.0/mu0_);
   }

   return coef;
}

// A spherical shell with constant permeability.  The sphere has inner
// and outer radii, center, and relative permeability specified on the
// command line and stored in ms_params_.
double magnetic_shell(const Vector &x)
{
   double r2 = 0.0;

   for (int i = 0; i < x.Size(); i++)
   {
      r2 += (x(i) - ms_params_(i))*(x(i) - ms_params_(i));
   }

   if ( sqrt(r2) >= ms_params_(x.Size()) &&
        sqrt(r2) <= ms_params_(x.Size()+1) )
   {
      return mu0_*ms_params_(x.Size()+2);
   }
   return mu0_;
}

// An annular ring of current density.  The ring has two axis end
// points, inner and outer radii, and a constant current in Amperes.
void current_ring(const Vector &x, Vector &j)
{

   // Current Density Function
   // Axis Start, Axis End, Inner Ring Radius,
   //                               Outer Ring Radius, and Total Current
   //                               of current ring (annulus)
   Vector cr_params_(5);
   //      mpirun -np 4 tesla -cr '0 0 -0.2 0 0 0.2 0.2 0.4 1'
   cr_params_(0) = 0.; 
   cr_params_(1) = 0.; 
   cr_params_(2) = 0.; 
   cr_params_(3) = 2.; 
   cr_params_(4) = 0.; 
   cr_params_(5) = 0.; 
   cr_params_(6) = 0.2; 
   cr_params_(7) = 0.4; 
   cr_params_(8) = 1.; 

   MFEM_ASSERT(x.Size() == 3, "current_ring source requires 3D space.");

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

// A Cylindrical Rod of constant magnetization.  The cylinder has two
// axis end points, a radius, and a constant magnetic field oriented
// along the axis.
void bar_magnet(const Vector &x, Vector &m)
{
   m.SetSize(x.Size());
   m = 0.0;

   Vector  a(x.Size());  // Normalized Axis vector
   Vector xu(x.Size());  // x vector relative to the axis end-point

   xu = x;

   for (int i=0; i<x.Size(); i++)
   {
      xu[i] -= bm_params_[i];
      a[i]   = bm_params_[x.Size()+i] - bm_params_[i];
   }

   double h = a.Norml2();

   if ( h == 0.0 )
   {
      return;
   }

   double  r = bm_params_[2*x.Size()];
   double xa = xu*a;

   if ( h > 0.0 )
   {
      xu.Add(-xa/(h*h),a);
   }

   double xp = xu.Norml2();

   if ( xa >= 0.0 && xa <= h*h && xp <= r )
   {
      m.Add(bm_params_[2*x.Size()+1]/h,a);
   }
}

// A Square Rod of rotating magnetized segments.  The rod is defined
// by a bounding box and a number of segments.  The magnetization in
// each segment is constant and follows a rotating pattern.
void halbach_array(const Vector &x, Vector &m)
{
   m.SetSize(x.Size());
   m = 0.0;

   // Check Bounding Box
   if ( x[0] < ha_params_[0] || x[0] > ha_params_[3] ||
        x[1] < ha_params_[1] || x[1] > ha_params_[4] ||
        x[2] < ha_params_[2] || x[2] > ha_params_[5] )
   {
      return;
   }

   int ai = (int)ha_params_[6];
   int ri = (int)ha_params_[7];
   int n  = (int)ha_params_[8];

   int i = (int)n * (x[ai] - ha_params_[ai]) /
           (ha_params_[ai+3] - ha_params_[ai]);

   m[(ri + 1 + (i % 2)) % 3] = pow(-1.0,i/2);
}

// To produce a uniform magnetic flux the vector potential can be set
// to ( By z, Bz x, Bx y).
void a_bc_uniform(const Vector & x, Vector & a)
{
   a.SetSize(3);
   a(0) = b_uniform_(1) * x(2);
   a(1) = b_uniform_(2) * x(0);
   a(2) = b_uniform_(0) * x(1);
}

// To produce a uniform magnetic field the scalar potential can be set
// to -z (or -y in 2D).
double phi_m_bc_uniform(const Vector &x)
{
   return -x(x.Size()-1);
}
