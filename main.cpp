#include <iostream>
#include <cmath>
#include <numbers>
#include <random>
#include <complex>

#include "vector.h"
#include "matrix.h"
#include "nonLinear.h"
#include "powerIterations.h"

double f(double x) {
    return (std::pow(x, 9) + std::numbers::pi)
        * cos(log(std::pow(x, 2) + 1)) /
        std::pow(std::numbers::e, std::pow(x, 2)) - x / 2022;
}

/*double f(double x) {
    return std::pow(2, x) - x - 3;
}*/

double df(double x) {
    return -2 * std::pow(std::numbers::e, -std::pow(x, 2)) * (std::pow(x, 9) + std::numbers::pi)
           * x * sin(log(std::pow(x, 2) + 1)) / (std::pow(x, 2) + 1)

           -2 * std::pow(std::numbers::e, -std::pow(x, 2)) * (std::pow(x, 9) + std::numbers::pi) * x
            * cos(log(std::pow(x, 2) + 1))

           +9 * std::pow(std::numbers::e, -std::pow(x, 2))
              * std::pow(x, 8) * cos(log(std::pow(x, 2) + 1)) - 1.0/2022;

}

std::vector<std::complex<double>> findR(Vector* lv) {
    //Vector l3(lv->length());
    std::vector<std::complex<double>> l3;
    for (int i = 0; i < lv->length(); i++) {
        double r = std::sqrt((lv[1](i) * lv[3](i)) - std::pow(lv[2](i), 2) / (lv[0](i) * lv[2](i) - std::pow(lv[1](i), 2)));
        double cos = (lv[3](i) + std::pow(r, 2) * lv[1](i)) / (2 * r * lv[2](i));
        double sin = 1 - std::pow(cos, 2);

        l3.emplace_back(r * cos, r * sin);
    }

    return l3;
}

Vector powerIteration(const Matrix& a, Vector u0, double h) {
    //Vector l1(u0.length());
    //Vector l2(u0.length());
    //Vector l3(u0.length());

    Vector v(u0.length());
    Vector u(u0.length());
    Vector um1(u.length());
    //prev(0) = -10;

    /*for (int i = 0; i < 100; i++) {
        v = (a * u0);
        u = v / v.norm();
        l1 = v / u0;
        u0 = u;
    }*/

    auto* lastVec = new Vector[4];
    auto* l1 = new Vector[4];
    auto* l2 = new Vector[4];
    std::vector<std::complex<double>> l3;
    //Vector l3(u0.length());
    //auto* l3 = new Vector[4];
    for (int i = 0; i < 100; i++) {
        u = (a * u0);
        //l1 = v / u0;
        switch (i%4) {
            case 0:
                lastVec[0] = u;
                break;
            case 1:
                lastVec[1] = u;
                break;
            case 2:
                lastVec[2] = u;
                break;
            case 3:
                lastVec[3] = u;

                l1[0] = lastVec[1] / lastVec[0];
                l1[1] = lastVec[2] / lastVec[1];
                l1[2] = lastVec[3] / lastVec[2];

                l2[0] = lastVec[2] / lastVec[0];
                l2[1] = lastVec[3] / lastVec[1];
                for(int i = 0; i < l2[0].length(); i++) {
                    l2[0](i) = std::sqrt(std::abs(l2[0](i)));
                    l2[1](i) = std::sqrt(std::abs(l2[1](i)));
                }

                l3 = findR(lastVec);


                if ((l2[1] - l2[0]).norm() < h) {
                    std::cout << "PICKED SQUARES METHOD" << std::endl;
                    std::cout << "EIGENVEC" << lastVec[3] + (l2[1] * lastVec[2]) << std::endl;
                    std::cout << "EIGENVAL" << l2[1] << std::endl;


                    return l2[1];
                } else if ((l1[2] - l1[1]).norm() < h) {
                    std::cout << "PICKED USUAL METHOD" << std::endl;
                    std::cout << "EIGENVEC" << lastVec[3] / lastVec[3].norm() << std::endl;
                    std::cout << "EIGENVAL" << l1[2] << std::endl;

                    return l1[2];

                }

                u = u / u.norm();
                break;
        }

        u0 = u;

    }

    /*std::cout << "----------USUAL----------" << std::endl;
    std::cout << "EIGENVEC" << lastVec[3] / lastVec[3].norm() << std::endl;
    std::cout << "EIGENVAL" << l1[2] << std::endl;

    std::cout << "----------SQUARE----------" << std::endl;
    std::cout << "EIGENVEC" << lastVec[3] + (l2[1] * lastVec[2]) << std::endl;
    std::cout << "EIGENVAL" << l2[1] << std::endl;*/

    /*std::cout << "----------COMPLEX----------" << std::endl;
    std::cout << "EIGENVEC" << lastVec[1] + lastVec[0] * l3 << std::endl;
    std::cout << "EIGENVAL" << l3.at(0).real() << std::endl;*/

    //std::cout << ((a * r) * r) / (r * r) << std::endl;

    return l1[2];
}

int main() {
    /*auto roots = NonLinear(f, -10, 10).solve();
    for (auto v : roots) {
        std::cout << v << std::endl;
    }*/

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
    //std::cout << m1 << std::endl;

    Matrix test{ {18, -8, -20}, {20, -10, -20}, {8, -8, -10} };
    //std::cout << m1 << std::endl;

    Vector r1(10, 1);
    Vector r2(20, 1);
    std::cout << r1 << std::endl;

    //std::cout << powerIteration(test, r1, 10e-8);
    std::cout << PowerIterations::powerIteration(m2, r2, 10e-8) << std::endl;

    return 0;
}
