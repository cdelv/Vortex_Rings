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
#include "diffusion.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "lambda2.h"
#include "utils.h"

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
  -
*/
#include "Integrals/connector.h"

// Macro for ring refinement
#define RADIUS (sqrt(sq(x) + sq(y) + sq(z)))

struct Config {
    // Finalization time.
    double tmax;

    // Box size.
    double L;

    // Ring Parameters: radius, Initial position in the z-direction,
    // Sigma of the gaussian (thickness of the ring), Magnitude of the vortex, and
    // Tolerance for the Laplace solver.
    // R=1 and Gamma=1 always.
    double R, Z0, a, Gamma, tol;

    // Fluid parameters
    double Re, viscosity;

    // AMR parameters
    int max_level, initial_level;
    double threshold;

    // File saving parameters
    char *path;
} conf;

// How often to save data.
double save_dt = 0.1;

// Initialize the values of the Config struct
void init_values(int argc, char *argv[]);

/*
  - curl: computes the curl of a vector field.
  - Integral: computes ....
*/
void curl(const vector v, vector curl);

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
p[back]   =  neumann(0.0);
p[top]    = dirichlet(0.0);
p[bottom] = dirichlet(0.0);
p[right]  = dirichlet(0.0);
p[left]   = dirichlet(0.0);

int main(int argc, char *argv[]) {
    init_values(argc, argv);
    init_grid ((int)pow(2, conf.initial_level));
    size(conf.L);
    X0 = Y0 = Z0 = -L0 / 2;
    const face vector muc[] = {conf.viscosity, conf.viscosity, conf.viscosity};
    mu = muc;
    run();
}

/*
  Initial Condition: we initialize the velocity from a vorticity field.
*/
event init(t = 0.0) {
    // Refine the ring
    refine (RADIUS < conf.R + 1.1 * conf.a && level < conf.max_level - 1);
    refine (RADIUS < conf.R + 4.5 * conf.a && level < conf.max_level - 2);
    unrefine(RADIUS > conf.R + 5.0 * conf.a && level > 4);
    unrefine(RADIUS > conf.R + 10.0 * conf.a && level > 3);
    boundary (all);

    if (pid() == 0)
        printf("\n %s ", "Computing Integral ... ");

    foreach () {
        u.x[] = U_x0(x,y,z,conf.Z0,conf.a);
        u.y[] = U_y0(x,y,z,conf.Z0,conf.a);
        u.z[] = U_z0(x,y,z,conf.Z0,conf.a);
        printf("%6f %6f %6f %g %g %g %d \n", x, y, z, u.x[], u.y[], u.z[], pid());

    }

    if (pid() == 0)
        printf("%s\n", "Done.");
}


/*
  Perform the AMR.
  - We use the velocity as the criteria for refinement.
*/
event adapt (i++) {
    astats s = adapt_wavelet ((scalar*) {u}, (double[]) {conf.threshold , conf.threshold, conf.threshold}, conf.max_level, 3);
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
  CHECK FOR FILE EXISTANCE
*/
event snapshots (t += save_dt) {
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
            W_Z.x[] = (z - conf.Z0) * W_vec.x[];
            W_Z2.x[] = pow(z - conf.Z0, 2) * W_vec.x[];
        }
    }

    // Paraview output
    char prefix[80];
    sprintf(prefix, "data_%03d_%06d", (int) t, i);
    output_htg((scalar *) {W_mag, PID, LEVEL}, (vector *) {u, W_vec, W_r, W_r2, W_Z, W_Z2}, conf.path, prefix, i, t);
}

event stop (t = conf.tmax);

void init_values(int argc, char *argv[]) {
    char *ptr;
    conf.tmax = strtod(argv[1], &ptr);
    conf.L = strtod(argv[2], &ptr);
    conf.R = 1.0;
    conf.Z0 = strtod(argv[3], &ptr);
    conf.a = strtod(argv[4], &ptr);
    conf.Gamma = 1.0;
    conf.tol = strtod(argv[5], &ptr);
    conf.Re = strtod(argv[6], &ptr);
    conf.viscosity = 1.0 / conf.Re;
    conf.initial_level = strtod(argv[7], &ptr);
    conf.max_level = strtod(argv[8], &ptr);
    conf.threshold = strtod(argv[9], &ptr);
    save_dt = strtod(argv[10], &ptr);
    conf.path = argv[11];
    if (pid() == 0) {
        printf("tmax = %g \n", conf.tmax);
        printf("L = %g \n", conf.L);
        printf("Z0 = %g \n", conf.Z0);
        printf("a = %g \n", conf.a);
        printf("tol = %g \n", conf.tol);
        printf("Re = %g \n", conf.Re);
        printf("initial_level = %d \n", conf.initial_level);
        printf("max_level = %d \n", conf.max_level);
        printf("threshold = %g \n", conf.threshold);
        printf("save_dt = %g \n", save_dt);
        printf("path = %s \n", conf.path);
    }
}

void Integral(double x, double y, double z, double vals[3]) {
    vals[0] = 0.0;
    vals[1] = 0.0;
    vals[2] = 0.0;

    // Create the command to be executed
    char command[256];
    sprintf (command, "./integral %g %g %g %g %g", conf.Z0, conf.a, x, y, z);
    printf("%s\n", command);

    // Execute the command
    FILE *cmd = popen(command, "r");

    // Get the command output
    char buf[256] = {0x0};
    while (fgets(buf, sizeof(buf), cmd) != NULL) {
        // Output is separated with ,
        // separate the string
        char *token = strtok(buf, ",");
        for (int i = 0; i < 3; ++i) {
            // store integral value
            vals[i] = atof(token);
            token = strtok(NULL, ",");
        }
    }

    // Close file
    pclose(cmd);
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