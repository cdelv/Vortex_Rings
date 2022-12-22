#include <cmath>
#include <iostream>
#include "vector.h"

const int points = 30;


#define Gauss

#ifdef Gauss
    #include "boost/math/quadrature/gauss_kronrod.hpp"
#else
    #include <boost/math/quadrature/sinh_sinh.hpp>
#endif

using namespace boost::math::quadrature;


/* The center of the ring is located in
    x = 0.0
    y = 0.0
    z = Z0

The ring radius R, by defect in the Basilisk program
    R = 1.0

The ring vorticity magnitude Gamma, by defect in the Basilisk program
    Gamma = 1.0

The ring thickness a
*/
struct Config{
    double R = 1.0;
    double Gamma = 1.0;
    double Z0;
    double a;
    double eps = std::sqrt(std::numeric_limits<double>::epsilon());
} conf;

double Integral_x(vector3D r);
double Integral_y(vector3D r);
double Integral_z(vector3D r);
void Vorticity(vector3D &x, vector3D &u);

int main(int argc, char *argv[])
{   
    conf.Z0 = atof(argv[1]);
    conf.a = atof(argv[2]);

    // Points to evaluate the integral (r)
    double x = atof(argv[3]);
    double y = atof(argv[4]);
    double z = atof(argv[5]);
    
    // W the vorticity and u the velocity
    // W = ∇ x u
    // ∇⋅u = 0
    // u = 1/4π ∫( W(r')x(r - r')/|r-r'|³ )d³r' 
    // prints u.x
    double u_x = Integral_x(vector3D(x,y,z));
    double u_y = Integral_y(vector3D(x,y,z));
    double u_z = Integral_z(vector3D(x,y,z));

    std::cout << std::fixed << std::setprecision(16) << u_x << "," << u_y << "," << u_z << ",";

    return 0;
}

 
void Vorticity(vector3D &x, vector3D &u)
{ 
    double theta = std::atan2(x.z()-conf.Z0, x.y());
    u.load(0.0,conf.R*std::cos(theta),conf.Z0 + conf.R*std::sin(theta));   
    double r = norm(u-x);                        
 
    if (r <= conf.a){
        u.load(0.0,-std::sin(theta),std::cos(theta));
        u *= (1-r/conf.a)*conf.Gamma;   
    } else
        u.load(0.0,0.0,0.0);
}

double Integral_x(vector3D r){

    #ifdef Gauss
        double x1 = -std::numeric_limits<double>::infinity();
        double x2 = std::numeric_limits<double>::infinity();
        double y1 = -std::numeric_limits<double>::infinity();
        double y2 = std::numeric_limits<double>::infinity();
        double z1 = -std::numeric_limits<double>::infinity();
        double z2 = std::numeric_limits<double>::infinity();
    #else
        sinh_sinh<double> integrator(4);
        double error;
        double L1;
    #endif

    vector3D cross, W, rr;

    auto f2 = [&](double x, double y, double z){
        //W(r')
        rr.load(x,y,z);
        Vorticity(rr,W);

        //r-r'
        rr = r - rr;

        //W(r') x (r-r')
        cross = W^rr;

        double norm = rr.norm2();
        // cross(coord)
        return cross.x()/std::pow(norm+conf.eps,3);
    };

    auto f1 = [&](double x, double y) { 
        auto g = [&](double z) { return f2(x,y,z); };
        #ifdef Gauss
            return gauss<double, points>::integrate(g, z1, z2);
        #else
            return integrator.integrate(g, error, &L1);
        #endif
    };

    auto f = [&](double x) { 
        auto g = [&](double y) { return f1(x, y); };
        #ifdef Gauss
            return gauss<double, points>::integrate(g, y1, y2);
        #else
            return integrator.integrate(g, error, &L1);
        #endif
    };
    #ifdef Gauss
        return 0.25*M_1_PI*gauss<double, points>::integrate(f, x1, x2);
    #else
        return 0.25*M_1_PI*integrator.integrate(f, error, &L1);
    #endif
}

double Integral_y(vector3D r){

    #ifdef Gauss
        double x1 = -std::numeric_limits<double>::infinity();
        double x2 = std::numeric_limits<double>::infinity();
        double y1 = -std::numeric_limits<double>::infinity();
        double y2 = std::numeric_limits<double>::infinity();
        double z1 = -std::numeric_limits<double>::infinity();
        double z2 = std::numeric_limits<double>::infinity();
    #else
        sinh_sinh<double> integrator(4);
        double error;
        double L1;
    #endif

    vector3D cross, W, rr;

    auto f2 = [&](double x, double y, double z){
        //W(r')
        rr.load(x,y,z);
        Vorticity(rr,W);

        //r-r'
        rr = r - rr;

        //W(r') x (r-r')
        cross = W^rr;

        double norm = rr.norm2();
        // cross(coord)
        return cross.y()/std::pow(norm+conf.eps,3);
    };

    auto f1 = [&](double x, double y) { 
        auto g = [&](double z) { return f2(x,y,z); };
        #ifdef Gauss
            return gauss<double, points>::integrate(g, z1, z2);
        #else
            return integrator.integrate(g, error, &L1);
        #endif
    };

    auto f = [&](double x) { 
        auto g = [&](double y) { return f1(x, y); };
        #ifdef Gauss
            return gauss<double, points>::integrate(g, y1, y2);
        #else
            return integrator.integrate(g, error, &L1);
        #endif
    };
    #ifdef Gauss
        return 0.25*M_1_PI*gauss<double, points>::integrate(f, x1, x2);
    #else
        return 0.25*M_1_PI*integrator.integrate(f, error, &L1);
    #endif
}

double Integral_z(vector3D r){

    #ifdef Gauss
        double x1 = -std::numeric_limits<double>::infinity();
        double x2 = std::numeric_limits<double>::infinity();
        double y1 = -std::numeric_limits<double>::infinity();
        double y2 = std::numeric_limits<double>::infinity();
        double z1 = -std::numeric_limits<double>::infinity();
        double z2 = std::numeric_limits<double>::infinity();
    #else
        sinh_sinh<double> integrator(4);
        double error;
        double L1;
    #endif

    vector3D cross, W, rr;

    auto f2 = [&](double x, double y, double z){
        //W(r')
        rr.load(x,y,z);
        Vorticity(rr,W);

        //r-r'
        rr = r - rr;

        //W(r') x (r-r')
        cross = W^rr;

        double norm = rr.norm2();
        // cross(coord)
        return cross.z()/std::pow(norm+conf.eps,3);
    };

    auto f1 = [&](double x, double y) { 
        auto g = [&](double z) { return f2(x,y,z); };
        #ifdef Gauss
            return gauss<double, points>::integrate(g, z1, z2);
        #else
            return integrator.integrate(g, error, &L1);
        #endif
    };

    auto f = [&](double x) { 
        auto g = [&](double y) { return f1(x, y); };
        #ifdef Gauss
            return gauss<double, points>::integrate(g, y1, y2);
        #else
            return integrator.integrate(g, error, &L1);
        #endif
    };
    #ifdef Gauss
        return 0.25*M_1_PI*gauss<double, points>::integrate(f, x1, x2);
    #else
        return 0.25*M_1_PI*integrator.integrate(f, error, &L1);
    #endif
}