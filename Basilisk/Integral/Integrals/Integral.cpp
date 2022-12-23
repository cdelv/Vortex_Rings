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

const int points = 100;


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
struct Config {
    double R = 1.0;
    double Gamma = 1.0;
    double Z0;
    double a;
    double eps = 1.0e-9;
} Conf;

void Vorticity(vector3D &x, vector3D &u)
{
    double theta = std::atan2(x.z() - Conf.Z0, x.y());
    u.load(0.0, Conf.R * std::cos(theta), Conf.Z0 + Conf.R * std::sin(theta)); // linear vortex
    double r = norm(u - x);

    // lineal ring
    /*if (r <= Conf.a) {
        u.load(0.0, -std::sin(theta), std::cos(theta)); // tangent vortex
        u *= (1 - r / Conf.a) * Conf.Gamma;
    } else
        u.load(0.0, 0.0, 0.0);*/

    //gausian ring
    u.load(0.0, -std::sin(theta), std::cos(theta)); // tangent vortex
    u *= (Conf.Gamma/(M_PI*std::pow(Conf.a,2)))*std::exp(-std::pow(r/Conf.a,2));
}

double Integral_x(double xx, double yy, double zz, double Z0, double a) {
    Conf.a = a;
    Conf.Z0 = Z0;
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
    return 0.25 * M_1_PI * gauss<double, points>::integrate(f, x1, x2);
#else
    return 0.25 * M_1_PI * integrator.integrate(f, error, &L1);
#endif
}

double Integral_y(double xx, double yy, double zz, double Z0, double a) {
    Conf.a = a;
    Conf.Z0 = Z0;
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
    return 0.25 * M_1_PI * gauss<double, points>::integrate(f, x1, x2);
#else
    return 0.25 * M_1_PI * integrator.integrate(f, error, &L1);
#endif
}

double Integral_z(double xx, double yy, double zz, double Z0, double a) {
    Conf.a = a;
    Conf.Z0 = Z0;
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
        return cross.z() * std::pow(norm + Conf.eps, -3);
    };

    auto f1 = [&](double x, double y) {
        auto g = [&](double z) { return f2(x, y, z); };
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
    return 0.25 * M_1_PI * gauss<double, points>::integrate(f, x1, x2);
#else
    return 0.25 * M_1_PI * integrator.integrate(f, error, &L1);
#endif
}