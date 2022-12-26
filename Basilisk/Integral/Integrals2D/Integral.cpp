#include <cmath>
#include "boost/math/quadrature/gauss_kronrod.hpp"

const int points = 15;
const int depth = 6;

double Ux = 0, Uy = 0, Uz = 0;

struct Config {
    double R = 1.0;
    double Gamma = 1.0;
    double Z0;
    double a;
    double eps = 1.0e-6;
    int n = 20;
} Conf;

double w(double rho, double z) {
    return Conf.Gamma * M_1_PI * std::pow(Conf.a, -2) * std::exp(-(std::pow((rho - Conf.R) / Conf.a, 2) + std::pow((z - Conf.Z0) / Conf.a, 2)) );
}

double dE_dk(double x) {
    if (x < 0.01) {
        double sum = 0;
        for (int j = 0; j < Conf.n; j++)
            sum += std::pow(std::legendre(2 * j + 2, 0.0), 2) * (2 * j + 2) * std::pow(x, 2 * j) / (2 * j + 1);
        return -M_PI_2 * sum;
    }
    else
        return (std::comp_ellint_2(x) - std::comp_ellint_1(x)) / std::pow(x, 2);
}

double G_rho(double rho, double z, double rhop, double zp) {
    double r1_2 = std::pow(rho - rhop, 2) + std::pow(z - zp, 2) + Conf.eps;
    double r2_2 = std::pow(rho + rhop, 2) + std::pow(z - zp, 2) + Conf.eps;
    double arg = 2.0 * std::sqrt(rho * rhop / r2_2);
    double E = std::comp_ellint_2(arg);
    double F = dE_dk(arg);
    return rhop * (z - zp) * M_1_PI * (E / r1_2 + 2 * F / r2_2) / std::sqrt(r2_2);
}

double G_z(double rho, double z, double rhop, double zp) {
    double r1_2 = std::pow(rho - rhop, 2) + std::pow(z - zp, 2) + Conf.eps;
    double r2_2 = std::pow(rho + rhop, 2) + std::pow(z - zp, 2) + Conf.eps;
    double arg = 2.0 * std::sqrt(rho * rhop / r2_2);
    double E = std::comp_ellint_2(arg);
    double F = dE_dk(arg);
    return rhop * M_1_PI * ((rhop - rho) * E / r1_2 - 2 * rho * F / r2_2) / std::sqrt(r2_2);
}

double IG_rho(double rho, double z) {
    double rho1 = 0.0;
    double rho2 = std::numeric_limits<double>::infinity();
    double z1 = -std::numeric_limits<double>::infinity();
    double z2 = std::numeric_limits<double>::infinity();
    double error, L1;

    auto Integrand = [&rho, &z](double rhop, double zp) {return w(rhop, zp) * G_rho(rho, z, rhop, zp);};

    auto f = [&](double rhop) {
        auto g = [&](double zp) {
            return Integrand(rhop, zp);
        };
        return boost::math::quadrature::gauss_kronrod<double, points>::integrate(g, z1, z2, depth, Conf.eps);
    };

    return boost::math::quadrature::gauss_kronrod<double, points>::integrate(f, rho1, rho2, depth, Conf.eps);
}

double IG_z(double rho, double z) {
    double rho1 = 0.0;
    double rho2 = std::numeric_limits<double>::infinity();
    double z1 = -std::numeric_limits<double>::infinity();
    double z2 = std::numeric_limits<double>::infinity();

    double error, L1;

    auto Integrand = [&rho, &z](double rhop, double zp) {return w(rhop, zp) * G_z(rho, z, rhop, zp);};

    auto f = [&](double rhop) {
        auto g = [&](double zp) {
            return Integrand(rhop, zp);
        };
        return boost::math::quadrature::gauss_kronrod<double, points>::integrate(g, z1, z2, depth, Conf.eps);
    };

    return boost::math::quadrature::gauss_kronrod<double, points>::integrate(f, rho1, rho2, depth, Conf.eps);
}

void compute_V0(double x, double y, double z, double a, double Z0) {
    Conf.a = a;
    Conf.Z0 = Z0;
    double rho = std::hypot(x, y);
    double theta = std::atan2(y, x);

    Uz = IG_z(rho, z);

    double V_rho = IG_rho(rho, z);

    Ux = V_rho * std::cos(theta);
    Uy = V_rho * std::sin(theta);
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

/*int main(int argc, char *argv[])
{
    double r, dr = 0.01;
    int n = 100;

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n; ++i)
    {
        IG_z(1.0 - 1.12 * Conf.a, 0.0);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Time : " << time_taken / n << " ms" << std::endl;


    std::cout << "a = " << Conf.a << std::endl;
    std::cout << 0.64 * M_1_PI * 0.5 / Conf.a << std::endl;
    std::cout << IG_z(1.0 - 1.12 * Conf.a, 0.0) << std::endl;

    for(r=0; r<2; r+=dr)
        std::cout << r << "\t" << IG_z(r, 0.0) << std::endl;

    return 0;
}*/