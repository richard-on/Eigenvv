#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <random>
#include <complex>
#include <chrono>

#include "vector.h"
#include "matrix.h"
#include "nonLinear.h"
#include "eigen.h"

double f(double x) {
    return (std::pow(x, 9) + M_PI)
        * cos(log(std::pow(x, 2) + 1)) /
        std::pow(M_E, std::pow(x, 2)) - x / 2022;
}

double df(double x) {
    return -2 * std::pow(M_E, -std::pow(x, 2)) * (std::pow(x, 9) + M_PI)
           * x * sin(log(std::pow(x, 2) + 1)) / (std::pow(x, 2) + 1)

           -2 * std::pow(M_E, -std::pow(x, 2)) * (std::pow(x, 9) + M_PI) * x
            * cos(log(std::pow(x, 2) + 1))

           +9 * std::pow(M_E, -std::pow(x, 2))
              * std::pow(x, 8) * cos(log(std::pow(x, 2) + 1)) - 1.0/2022;

}

void nonLinearExample() {

    NonLinear nl(f, -10, 10, 0.1);

    auto roots = nl.getRootIntervals();
    for (auto v : roots) {
        Solution bisected = nl.bisect(v, 1e-4);
        Solution nl2 = nl.newton(v);

        std::cout << "Narrowed interval after " << bisected.iterations << " bisect iterations" << std::endl;
        std::cout << "Using interval [" << bisected.interval.start << "; " << bisected.interval.end << "]" << std::endl;
        std::cout << "Found root after " << nl2.iterations << " newton iterations" << std::endl;
        std::cout << "Root: " << nl2.root << std::endl << std::endl;
    }

}

int main() {
    nonLinearExample();

    Matrix m1{{1,-2,1,0,-1,1,-2,2,0,-2},
              {0,2,0,0,2,1,-1,-1,-1,-2},
              {0,1,0,-1,1,-1,0,-1,1,-1},
              {-2,-1,2,-1,0,0,0,0,1,0},
              {1,-2,0,1,0,-2,-1,0,2,2},
              {-2,-2,0,-2,0,1,1,-2,1,1},
              {-1,-2,-1,-1,-2,-1,-2,1,-1,2},
              {-2,1,2,-2,0,2,1,-1,-2,2},
              {0,1,0,1,1,-2,2,0,1,1},
              {0,0,2,-1,-1,0,-2,2,-1,-1}};

    Matrix m2{{ -1, 1, -1, 0, -1, 0, -1, 1, 1, -1, 0, -1, -1, 1, 0, 0, 1, 1, 1, 1 },
            { -1, 0, -1, 1, -1, 0, 0, 0, 0, -1, 0, 0, -1, 1, 0, -1, 1, -1, -1, 0 },
            { 1, 0, -1, 1, 0, 1, -1, -1, -1, 0, -1, -1, 1, -1, 1, 1, -1, 1, -1, 0 },
            { -1, 1, 0, 0, -1, 0, 0, -1, 0, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 0 },
            { 1, 0, -1, 0, 0, -1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, -1, 0, 0, 1 },
            { 0, 0, 0, 0, -1, 1, 1, 0, 0, 1, 1, 0, -1, 0, 1, 1, 0, 1, 0, 0 },
            { -1, 0, 1, 1, 1, -1, -1, 0, -1, 1, -1, -1, -1, 0, -1, 0, 0, 0, -1, 1 },
            { 0, 0, -1, -1, 0, 1, 1, 1, 1, -1, 0, 0, -1, 1, 1, 1, 1, 0, 0, -1 },
            { 0, 0, 1, 1, 0, 1, 1, 0, 1, -1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1 },
            { 0, -1, 0, 0, 1, 0, -1, 0, -1, 0, -1, 0, -1, 0, 1, -1, 0, 0, 1, 1 },
            { 1, -1, 1, -1, -1, -1, 1, 0, -1, 0, 1, 1, -1, 0, 1, 1, 1, 0, 0, 0 },
            { 0, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 1, -1, -1, 0, -1, 1, 1, -1 },
            { -1, -1, -1, -1, 0, 1, -1, 0, 0, -1, 0, 0, 0, 1, 1, 0, 0, 0, -1, 0 },
            { -1, 0, 1, 0, -1, 0, 0, 1, -1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 0 },
            { 1, -1, 0, -1, -1, 0, -1, -1, 0, 0, 1, 0, 1, 1, -1, 1, 0, 0, -1, 0 },
            { -1, -1, 1, 0, -1, 1, 1, -1, 1, 0, 0, -1, 1, -1, -1, 0, 0, 1, 1, 1 },
            { 0, 0, -1, 0, 0, 0, 0, -1, 1, 1, 0, -1, 1, -1, 0, 0, 0, -1, -1, 1 },
            { -1, 0, -1, -1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, 0, -1, 0, -1 },
            { -1, 0, 1, 0, 0, 0, 0, -1, 1, -1, 1, -1, 0, -1, -1, 1, 0, 1, 0, 0 },
            { 0, -1, -1, 1, -1, 1, -1, -1, -1, 1, 1, -1, 0, -1, -1, 0, 1, 0, -1, -1 } };

    std::cout << eigen::powerIteration(m1) << std::endl;
    std::cout << eigen::powerIteration(m2) << std::endl;

    std::cout << eigen::QR(m1) << std::endl;
    std::cout << eigen::QR(m2) << std::endl;


    return 0;
}
