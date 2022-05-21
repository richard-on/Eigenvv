#include <iostream>
#include <cmath>
#include <numbers>
#include <random>

#include "vector.h"
#include "matrix.h"
#include "nonLinear.h"

double f(double x) {
    return (std::pow(x, 9) + std::numbers::pi)
        * cos(log(std::pow(x, 2) + 1)) /
        std::pow(std::numbers::e, std::pow(x, 2)) - x / 2022;
}

double df(double x) {
    return -2 * std::pow(std::numbers::e, -std::pow(x, 2)) * (std::pow(x, 9) + std::numbers::pi)
           * x * sin(log(std::pow(x, 2) + 1)) / (std::pow(x, 2) + 1)

           -2 * std::pow(std::numbers::e, -std::pow(x, 2)) * (std::pow(x, 9) + std::numbers::pi) * x
            * cos(log(std::pow(x, 2) + 1))

           +9 * std::pow(std::numbers::e, -std::pow(x, 2))
              * std::pow(x, 8) * cos(log(std::pow(x, 2) + 1)) - 1.0/2022;

}

Vector powerIteration(const Matrix& a, Vector r0, double h) {
    Vector r(10);
    Vector prev(10);
    //prev(0) = -10;

    for (int i = 0; i < 100; i++) {
        r = (a * r0) / (a * r0).norm();
        r0 = r;
    }

    /*while (std::abs((prev - r)(0)) > h) {
        r = (a * r0) / (a * r0).norm();
        prev = r0;
        r0 = r;
    }*/
    std::cout << (a * r).norm() << std::endl;

    std::cout << ((a * r) * r) / (r * r) << std::endl;

    return r;
}

int main() {
    for (auto v : NonLinear(f, -10, 10).solve()) {
        std::cout << v << std::endl;
    }

    /*std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(-2, 5);
    Matrix a = Matrix(10, data);
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            a(i, j) = distr(gen);
        }
    }
    std::cout << a << std::endl;*/


    /*Vector r1{1,-1,1,1,1,1,1,1,1,1};
    std::cout << r1 << std::endl;

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
    std::cout << m1 << std::endl;*/


    //std::cout << powerIteration(a, r0, 10e-8);

    return 0;
}
