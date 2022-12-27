#include "Integral.h"

/*
g++ -fpic -shared test.cpp -o libtest.so
g++ -fpic -shared connector.cpp -L. -ltest -o libconnector.so
export LD_LIBRARY_PATH=`pwd`
qcc main.c -L. -lconnector -lstdc++ -o c_aaa -lm
*/

#include <cmath>
#include <iostream>
#include "vector.h"

const int points = 15;
const int depth = 8;

#define Gauss

#ifdef Gauss
#include "boost/math/quadrature/gauss_kronrod.hpp"
#else
#include <boost/math/quadrature/sinh_sinh.hpp>
#endif

using namespace boost::math::quadrature;

double Ux = 0, Uy = 0, Uz = 0;

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
struct Config {
    double R = 1.0;
    double Gamma = 1.0;
    double Z0;
    double a;
    double eps = 1.0e-6;
} Conf;

double w(double rho, double z) {
    return Conf.Gamma * M_1_PI * std::pow(Conf.a, -2) * std::exp(-(std::pow((rho - Conf.R) / Conf.a, 2) + std::pow((z - Conf.Z0) / Conf.a, 2)) );
}

void Vorticity(vector3D &x, vector3D &u)
{
    double theta = std::atan2(x.y(), x.x());
    double rho = std::hypot(x.x(), x.y()) - Conf.R;
    double z = x.z() - Conf.Z0;

    //gausian ring
    u.load(-std::sin(theta), std::cos(theta), 0.0); // tangent vortex
    u *= Conf.Gamma * M_1_PI * std::pow(Conf.a, -2) * std::exp(-std::pow(rho / Conf.a, 2) - std::pow(z / Conf.a, 2));
}

double Integral_x(double xx, double yy, double zz) {
    vector3D r(xx, yy, zz);
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

    vector3D cross(0.0, 0.0, 0.0), W(0.0, 0.0, 0.0), rr(0.0, 0.0, 0.0), rrr(0.0, 0.0, 0.0);

    auto f2 = [&](double x, double y, double z) {
        //W(r')
        rrr.load(x, y, z);
        Vorticity(rrr, W);

        //r-r'
        rr = r - rrr;

        //W(r') x (r-r')
        cross = W ^ rr;

        double norm = rr.norm();

        // cross(coord)
        return cross.x() * std::pow(norm + Conf.eps, -3);
    };

    auto f1 = [&](double x, double y) {
        auto g = [&](double z) { return f2(x, y, z); };
#ifdef Gauss
        return gauss_kronrod<double, points>::integrate(g, z1, z2, depth, Conf.eps);
#else
        return integrator.integrate(g, error, &L1);
#endif
    };

    auto f = [&](double x) {
        auto g = [&](double y) { return f1(x, y); };
#ifdef Gauss
        return gauss_kronrod<double, points>::integrate(g, y1, y2, depth, Conf.eps);
#else
        return integrator.integrate(g, error, &L1);
#endif
    };
#ifdef Gauss
    return 0.25 * M_1_PI * gauss_kronrod<double, points>::integrate(f, x1, x2, depth, Conf.eps);
#else
    return 0.25 * M_1_PI * integrator.integrate(f, error, &L1);
#endif
}

double Integral_y(double xx, double yy, double zz) {
    vector3D r(xx, yy, zz);
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

    vector3D cross(0.0, 0.0, 0.0), W(0.0, 0.0, 0.0), rr(0.0, 0.0, 0.0), rrr(0.0, 0.0, 0.0);

    auto f2 = [&](double x, double y, double z) {
        //W(r')
        rrr.load(x, y, z);
        Vorticity(rrr, W);

        //r-r'
        rr = r - rrr;

        //W(r') x (r-r')
        cross = W ^ rr;

        double norm = rr.norm();

        // cross(coord)
        return cross.y() * std::pow(norm + Conf.eps, -3);
    };

    auto f1 = [&](double x, double y) {
        auto g = [&](double z) { return f2(x, y, z); };
#ifdef Gauss
        return gauss_kronrod<double, points>::integrate(g, z1, z2, depth, Conf.eps);
#else
        return integrator.integrate(g, error, &L1);
#endif
    };

    auto f = [&](double x) {
        auto g = [&](double y) { return f1(x, y); };
#ifdef Gauss
        return gauss_kronrod<double, points>::integrate(g, y1, y2, depth, Conf.eps);
#else
        return integrator.integrate(g, error, &L1);
#endif
    };
#ifdef Gauss
    return 0.25 * M_1_PI * gauss_kronrod<double, points>::integrate(f, x1, x2, depth, Conf.eps);
#else
    return 0.25 * M_1_PI * integrator.integrate(f, error, &L1);
#endif
}

double Integral_z(double xx, double yy, double zz) {
    vector3D r(xx, yy, zz);
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

    vector3D cross(0.0, 0.0, 0.0), W(0.0, 0.0, 0.0), rr(0.0, 0.0, 0.0), rrr(0.0, 0.0, 0.0);

    auto f2 = [&](double x, double y, double z) {
        //W(r')
        rrr.load(x, y, z);
        Vorticity(rrr, W);

        //r-r'
        rr = r - rrr;

        //W(r') x (r-r')
        cross = W ^ rr;

        double norm = rr.norm();

        // W(r') x (r-r')/|r-r'|^3
        return cross.z() * std::pow(norm + Conf.eps, -3);
    };

    auto f1 = [&](double x, double y) {
        auto g = [&](double z) { return f2(x, y, z); };
#ifdef Gauss
        return gauss_kronrod<double, points>::integrate(g, z1, z2, depth, Conf.eps);
#else
        return integrator.integrate(g, error, &L1);
#endif
    };

    auto f = [&](double x) {
        auto g = [&](double y) { return f1(x, y); };
#ifdef Gauss
        return gauss_kronrod<double, points>::integrate(g, y1, y2, depth, Conf.eps);
#else
        return integrator.integrate(g, error, &L1);
#endif
    };
#ifdef Gauss
    return 0.25 * M_1_PI * gauss_kronrod<double, points>::integrate(f, x1, x2, depth, Conf.eps);
#else
    return 0.25 * M_1_PI * integrator.integrate(f, error, &L1);
#endif
}

void compute_V0(double x, double y, double z, double a, double Z0)
{
    Conf.a = a;
    Conf.Z0 = Z0;

    Ux = Integral_x(x, y, z);
    Uy = Integral_y(x, y, z);
    Uz = Integral_z(x, y, z);
}

double get_Vx(void) {
    return Ux;
}

double get_Vy(void) {
    return Uy;
}

double get_Vz(void) {
    return Uz;
}