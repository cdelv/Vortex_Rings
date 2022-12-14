//#include "grid/multigrid3D.h"
#include "grid/octree.h"
#include "poisson.h"
#include "navier-stokes/centered.h"
#include "output_htg.h"
#include "lambda2.h"

double Gamma = 1.0;
double aa = 0.2;
double R = 1.0;
double L = 0.0;

double W(double r, double z){
	double scale = Gamma/(M_PI*pow(aa,2));
	double T1 = pow(r-R,2);
	double T2 = pow(z-L,2);
	return scale*exp(-(T1+T2)/pow(aa,2));
}

double curl_w0_x(double x, double y, double z){
	double r = hypot(x,y)+1.0e-8;
	double T1 = 2.0*(z-L)/(r*pow(aa,2));
	return W(r,z)*T1*x;
}

double curl_w0_y(double x, double y, double z){
	double r = hypot(x,y)+1.0e-8;
	double T1 = 2.0*(z-L)/(r*pow(aa,2));
	return W(r,z)*T1*y;
}

double curl_w0_z(double x, double y, double z){
	double r = hypot(x,y)+1.0e-8;
	double T1 = 1.0/r;
	double T2 = z*(r-R)/pow(aa,2);
	return W(r,z)*(T1-T2);
}

int maxlevel = 8;

int main() {
  size (16.0);
  X0 = Y0 = Z0 = -L0/2;
  init_grid (64);

  scalar vx0[], vy0[], vz0[], bx[], by[], bz[];

  vx0[left] = dirichlet(0.0);
  vy0[left] = dirichlet(0.0);
  vz0[left] = dirichlet(0.0);

  foreach(){
  	bx[] = curl_w0_x(x,y,z);
  	by[] = curl_w0_y(x,y,z);
  	bz[] = curl_w0_z(x,y,z);
  }

  struct Poisson p;
  p.a = vx0;
  p.b = bx;
  p.tolerance = 0.001;
  poisson(p);

  p.a = vy0;
  p.b = by;
  poisson(p);

  p.a = vz0;
  p.b = bz;
  poisson(p);

  foreach(){
  	u.x[] = vx0[];
  	u.y[] = vy0[];
  	u.z[] = vz0[];
  }

  //output_ppm(a);
  // Paraview
  scalar l2[];
  lambda2 (u, l2); // vorticity.

  char path[]="htg"; // no slash at the end!!
  char prefix[20];
  sprintf(prefix, "data_%03d_%06d", (int) 0.1, 0);
  output_htg((scalar*){bx,by,bz,l2},(vector *){u}, path, prefix, 0, 0.1);
}



