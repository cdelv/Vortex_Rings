#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <iostream>

int main(int argc, char *argv[])
{
    using namespace boost::math::quadrature;

    double x1 = 0;
    double x2 = 2;
    double y1 = 0;
    double y2 = 3;
    double z1 = 0;
    double z2 = 4;


    auto f2 = [](double x, double y, double z){return 1/(x*x);};

    auto f1 = [&](double x, double y) { 

        auto g = [&](double z) { return f2(x,y,z); };

        return gauss_kronrod<double, 61>::integrate(g, z1, z2, 5);
    };

    auto f = [&](double x) { 

        auto g = [&](double y) { return f1(x, y); };

        return gauss_kronrod<double, 61>::integrate(g, y1, y2, 5);
    };

    double error;
    double Q = gauss_kronrod<double, 61>::integrate(f, x1, x2, 5, 1e-9, &error);
    std::cout << Q << ", error estimated at " << error <<std::endl;

    return 0;
}