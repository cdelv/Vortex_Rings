// Compile with
// qcc -source -D_MPI=1 -O2 two_rings.c
// mpicc -Wall -O2 -std=c99 _two_rings.c -lm -lmpi -L$BASILISK/gl -I$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa
// mpirun -np 4 ./a.out

/*
  This ones are basilisk includes.
  see http://basilisk.fr/Front%20Page
*/
#include "grid/octree.h"
#include "poisson.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "lambda2.h"
#include "utils.h"

/*
  To limit simulation time and to restart simulation if possible

#include "maxruntime.h"
# define MAXTIME 0.05*/

/*
  This include is for Paraview visualization. We took it from Sander Sandbox.
  Thank you very much for your help, Maximilian Sander.
  See http://basilisk.fr/sandbox/sander/output_htg.h.

  This exports htg (Hyper Tree Grid) data from the simulation to Paraview.
  We use this because Basilisk uses octrees for the AMR.

  Known HTG Problems
  - HyperTreeGridToDualGrid stops working when advancing the time step (Paraview related).
  - Contour Filter does not work on HTG with only one tree (like this exporter creates) (Paraview related).
  - x-z Axis is swapped (3D), x-y Axis is swapped (2D).
*/
#include "output_htg.h"

/*
  Time Variables:
  - tend: Finalization time.
*/
double tend = 20. + 0.1;

/*
  Fluid Parameters: CHANGE FOR VISCOSITY
  - Re: Reynolds number
*/
double Re = 1750.0;

/*
  Vortex Ring Parameters:
  -Gamma: Magnitude of the vortex ring.
  -aa: Sigma of the gaussian (thickness of the ring).
  -RR: Vortex radius.
  -LL: Initial position in the z-direction.
  -ttol: Tolerance for the Laplace solver.
*/
double Gamma = 1.0;
double aa = 0.2;
double RR = 1.5;
double LL = 10.0;
double ttol = 1.0e-8;

/*
  AMR variables:
  - maxlevel: maximum refinement depth in the tree.
  - ue: threshold for the refinement.
*/
int maxlevel = 8;
double ue = 0.008;

/*
  For vorticity initial conditions.
  - W: Scale factor of the initial vorticity curl.
  - curl_w0_x: X component of the initial vorticity curl.
  - curl_w0_y: Y component of the initial vorticity curl.
  - curl_w0_z: Z component of the initial vorticity curl.
*/
double W(double r, double z);
double curl_w0_x(double x, double y, double z);
double curl_w0_y(double x, double y, double z);
double curl_w0_z(double x, double y, double z);

/*
  Velocity boundary conditions:
  The order of the boundaries in 3D is
  { right,  left,   top, bottom, front, back  }
  { X max, X min, Y max,  Y min, Z max, Z min }
*/
u.n[right]  = neumann(0.0);
u.n[left]   = neumann(0.0);
u.n[top]    = neumann(0.0);
u.n[bottom] = neumann(0.0);
u.n[front]  = neumann(0.0);
u.n[back]   = neumann(0.0);

/*
  Pressure Boundary Conditions
*/
p[front]  = neumann(0.0);
p[back]  =  neumann(0.0);
p[top]    = dirichlet(0.0);
p[bottom] = dirichlet(0.0);
p[right]  = dirichlet(0.0);
p[left]   = dirichlet(0.0);

/*
  Laplace Boundary Conditions
  The Dirichlet condition is necessary for convergence.
  It's the farthest boundary from the ring.
*/
scalar vx0[], vy0[], vz0[];
scalar bx[], by[], bz[]; // forcing terms

vx0[back] = dirichlet(0.0);
vy0[back] = dirichlet(0.0);
vz0[back] = dirichlet(0.0);

/*
  Function for computing the curl of a vector field.
*/
void curl(const vector v, vector curl);

int main() {
    init_grid (64);
    size(32.0);
    X0 = Y0 = Z0 = -L0 / 2;
    const face vector muc[] = {1. / Re, 1. / Re, 1. / Re};
    mu = muc;
    run();
}

/*
  Initial Condition: we initialize the velocity from a vorticity field.
*/
event init(t = 0.0) {
    // Restart simulation from previus run
    //if (!restore (file = "restart")){
        // Initialize the forcing terms of the Poisson equation
        foreach () {
            bx[] = curl_w0_x(x, y, z);
            by[] = curl_w0_y(x, y, z);
            bz[] = curl_w0_z(x, y, z);
        }

        // Solve Poisson equation
        poisson (vx0, bx, tolerance = ttol);
        poisson (vy0, by, tolerance = ttol);
        poisson (vz0, bz, tolerance = ttol);

        // Initialize velocity
        foreach () {
            u.x[] = vx0[];
            u.y[] = vy0[];
            u.z[] = vz0[];
        }
    //}
}

/*
  Perform the AMR.
  - We use the velocity as the criteria for refinement.
*/
event adapt (i++) {
    astats s = adapt_wavelet ((scalar*) {u}, (double[]) {1.6 * ue, ue, ue}, maxlevel);
    fprintf (stderr, "# Time %3f step %d -> refined %d cells, coarsened %d cells\n", t, i, s.nf, s.nc);
}

/*
  Paraview Output:
  - W_vec: Vectorial vorticity field (W)
  - W_mag: Magnitude of the vorticity 
  - W_r: W*r where r = sqrt(x^2 + y^2)
  - W_r2: W*r^2 where r = sqrt(x^2 + y^2)
  - W_Z: W*(z-L) where L is the initial position of the ring
  - W_Z2: W*(z-L)^2 where L is the initial position of the ring
  - PID: Process ID of the CPU in charge of that cell
  - LEVEL: Level of refinement. The number of cells is 2^LEVEL
*/
event snapshots (t += 0.2) {
    vector W_vec[], W_r[], W_r2[], W_Z[], W_Z2[];
    scalar W_mag[], PID[], LEVEL[];

    curl(u, W_vec);

    foreach () {
        double r = hypot(x, y);
        W_mag[] = norm(W_vec);
        PID[] = pid();
        LEVEL[] = level;
        foreach_dimension() {
            W_r.x[] = r * W_vec.x[];
            W_r2.x[] = pow(r , 2) * W_vec.x[];
            W_Z.x[] = (z - LL) * W_vec.x[];
            W_Z2.x[] = pow(z - LL, 2) * W_vec.x[];
        }
    }

    // Paraview output
    char path[] = "htg"; // no slash at the end!!
    char prefix[80];
    sprintf(prefix, "data_%03d_%06d", (int) t, i);
    output_htg((scalar *) {W_mag, PID, LEVEL}, (vector *) {uf, W_vec, W_r, W_r2, W_Z, W_Z2}, path, prefix, i, t);
}

event stop (t = tend);

double W(double r, double z) {
    double scale = Gamma / (M_PI * pow(aa, 2));
    double T1 = pow(r - RR, 2);
    double T2 = pow(z - LL, 2);
    return scale * exp(-(T1 + T2) / pow(aa, 2));
}

double curl_w0_x(double x, double y, double z) {
    double r = hypot(x, y) + 1.0e-6;
    double T1 = 2.0 * (z - LL) / (r * pow(aa, 2));
    return W(r, z) * T1 * x;
}

double curl_w0_y(double x, double y, double z) {
    double r = hypot(x, y) + 1.0e-6;
    double T1 = 2.0 * (z - LL) / (r * pow(aa, 2));
    return W(r, z) * T1 * y;
}

double curl_w0_z(double x, double y, double z) {
    double r = hypot(x, y) + 1.0e-6;
    double T1 = 1.0 / r;
    double T2 = z * (r - RR) / pow(aa, 2);
    return W(r, z) * (T1 - T2);
}

void curl(const vector v, vector curl) {
    vector dvx[], dvy[], dvz[];
    scalar vx[], vy[], vz[];

    foreach () {
        vx[] = v.x[];
        vy[] = v.y[];
        vz[] = v.z[];
    }

    gradients ({vx}, {dvx});
    gradients ({vy}, {dvy});
    gradients ({vz}, {dvz});

    foreach () {
        curl.x[] = dvz.y[] - dvy.z[];
        curl.y[] = dvx.z[] - dvz.x[];
        curl.z[] = dvy.x[] - dvx.y[];
    }
}