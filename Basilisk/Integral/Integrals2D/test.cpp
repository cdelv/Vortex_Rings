#include "Integral.h"
#include <iostream>

int main(int argc, char *argv[])
{
    double x, y, z;
    double a = 0.01, Z0 = 0;
    compute_V0(0, 0, 0, a, Z0);


    std::cout << w(1,0) << "\t" << w(1+3*a,0) << "\n";

    
    return 0;
}