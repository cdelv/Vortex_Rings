// Compile with
// qcc -O2 -std=c99 -D_MPI=1 -I$(Paraview_Out) Circle_Ring.c -lm -lmpi -L$(pwd)/Integrals -lconnector -lstdc++ -Wl,-rpath,$(pwd)/Integrals
// export LD_LIBRARY_PATH=$(pwd)/Integrals
// mpirun -np 4 ./a.out

/*
  This ones are basilisk includes.
  see http://basilisk.fr/Front%20Page
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
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
  - Integrals3D/connector.h:
  - Integrals2D/connector.h:

  DONT FORGET TO CHENGE THIS IN CONFIG.MK AND LD_LIBRARY_PATHENV VARIABLE TOO
*/
//#include "Integrals3D/connector.h"
#include "Integrals2D/connector.h"

// Macro for ring refinement
#define R2 (sq(x) + sq(y))

struct Config {
    // Finalization time.
    double tmax;

    // Box size.
    double L;

    // Ring Parameters: radius, Initial position in the z-direction,
    // Sigma of the gaussian (thickness of the ring) and Magnitude of the vortex
    // R=1 and Gamma=1 always.
    double R, Z0, a, Gamma;

    // Fluid parameters
    double Re, viscosity;

    // AMR parameters
    double ns, threshold;
    int max_level, initial_level, min_level;
    double CFL;

    // File saving parameters
    char *path;
} conf;

// How often to save data.
double save_dt = 0.1;

/*
  - init_values: initialize the values of the Config struct and edits save_dt
  - curl: computes the curl of a vector field.
*/
void init_values(int argc, char *argv[]);
void curl(const vector v, vector curl);
double compute_threshold(double ns, double t);

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
    X0 = Y0 = -L0 / 2;
    Z0 = -conf.Z0;
    conf.Z0 = 0.0;
    const face vector muc[] = {conf.viscosity, conf.viscosity, conf.viscosity};
    mu = muc;
    run();
}

/*
  Initial Condition: we initialize the velocity from a vorticity field. We also refine the mesh
  with the vorticity before computing the integral.
*/
event init(t = 0.0) {
    // Set up the time step condition
    CFL = conf.CFL;
    // Compute Initial Vorticity.
    scalar W_mag[];
    compute_U0(0, 0, 0, conf.a, conf.Z0);
    foreach ()
        W_mag[] = W_0(hypot(x, y), z);

    // Refine the Ring Based on Vorticity.
    // As Vorticity is Gaussian, the Threshold is the Value at ns Standar Deviations.
    for (int ii = 0; ii < conf.initial_level; ++ii)
        adapt_wavelet ((scalar*) {W_mag}, (double[]) {W_0(conf.R + conf.ns * conf.a, conf.Z0)}, conf.max_level, conf.min_level + 1);

    printf("Number of elements in pid = %d: N = %d\n", pid(), grid->n);

    if (pid() == 0)
        printf("\n%s\n", "Computing The Integral ... ");

    // Compute Initial Velocity from the Vorticity.
    foreach () {
        compute_U0(x, y, z, conf.a, conf.Z0);
        u.x[] = U_x0();
        u.y[] = U_y0();
        u.z[] = U_z0();
        //printf("%6f %6f %6f %g %g %g %d \n", x, y, z, u.x[], u.y[], u.z[], pid());
    }

    if (pid() == 0)
        printf("%s\n", "Done.");
}

/*
  Perform the AMR.
  - We use the velocity as the criteria for refinement.
*/
event adapt (i++) {
    //conf.threshold = 0.5*exp(-conf.R*M_1_PI*t/(conf.a*conf.Re))*8.0e-4;
    conf.threshold = 2.0e-4/(1.0 + 4.0*M_PI*conf.a*t/(conf.Re*conf.R));

    scalar logu [];

    foreach()
        logu[] = norm(u);

    astats s = adapt_wavelet ((scalar*) {logu}, (double[]) {conf.threshold}, conf.max_level, conf.min_level);
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
  TO DO: CHECK FOR FILE EXISTANCE
*/
event snapshots (t += save_dt) {
    vector W_vec[];
    scalar W_mag[];

    curl(u, W_vec);

    foreach () {
        W_mag[] = norm(W_vec);
    }

    // Paraview output
    char prefix[80];
    sprintf(prefix, "data_%03d_%06d", (int) t, i);
    output_htg((scalar *) {W_mag}, (vector *) {u, W_vec}, conf.path, prefix, i, t);
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

    conf.Re = strtod(argv[5], &ptr);
    conf.viscosity = 1.0 / conf.Re;

    conf.ns = strtod(argv[6], &ptr);
    conf.initial_level = strtod(argv[7], &ptr);
    conf.max_level = strtod(argv[8], &ptr);
    conf.min_level = strtod(argv[9], &ptr);
    conf.CFL = strtod(argv[10], &ptr);

    save_dt = strtod(argv[11], &ptr);
    conf.path = argv[12];

    // print to comand line the configuration and
    // save a file with the configuration to conf.path.
    if (pid() == 0) {
        printf("tmax = %g \n", conf.tmax);
        printf("L = %g \n", conf.L);
        printf("R = %g \n", conf.R);
        printf("Z0 = %g \n", conf.Z0);
        printf("a = %g \n", conf.a);
        printf("Gamma = %g \n", conf.Gamma);
        printf("Re = %g \n", conf.Re);
        printf("ns = %g \n", conf.ns);
        printf("initial_level = %d \n", conf.initial_level);
        printf("max_level = %d \n", conf.max_level);
        printf("min_level = %d \n", conf.min_level);
        printf ("CFL = %g \n", conf.CFL);
        printf("save_dt = %g \n", save_dt);
        printf("path = %s \n", conf.path);
        char name[256];
        sprintf (name, "%s/conf.txt", conf.path);
        FILE * file;
        file = fopen (name, "w");
        fprintf (file, "tmax = %g \n", conf.tmax);
        fprintf (file, "L = %g \n", conf.L);
        fprintf (file, "R = %g \n", conf.R);
        fprintf (file, "Z0 = %g \n", conf.Z0);
        fprintf (file, "a = %g \n", conf.a);
        fprintf (file, "Gamma = %g \n", conf.Gamma);
        fprintf (file, "Re = %g \n", conf.Re);
        fprintf (file, "ns = %g \n", conf.ns);
        fprintf (file, "initial_level = %d \n", conf.initial_level);
        fprintf (file, "max_level = %d \n", conf.max_level);
        fprintf (file, "min_level = %d \n", conf.min_level);
        fprintf (file, "CFL = %g \n", conf.CFL);
        fprintf (file, "save_dt = %g \n", save_dt);
        fprintf (file, "path = %s \n", conf.path);
        time_t t = time(NULL);
        struct tm tm = *localtime(&t);
        fprintf(file, "%d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        fclose(file);
    }
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
