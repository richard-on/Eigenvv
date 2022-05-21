#include <complex>
#include <iostream>
#include "powerIterations.h"
#include <array>


Vector PowerIterations::powerIteration(const Matrix &a, Vector u0, double h, int maxIteration) {
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
    //std::array<Vector, 4> l1;
    //std::array<Vector, 4> l2;
    auto* l1 = new Vector[3];
    auto* l2 = new Vector[2];
    std::vector<std::complex<double>> l3(u0.length(), 0);
    std::vector<std::complex<double>> l3_old(u0.length(), 0);
    //Vector l3(u0.length());
    //auto* l3 = new Vector[4];
    for (int i = 0; i < maxIteration; i++) {
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

                l3 = findComplexL(lastVec);

                double delta_prev = 0;
                double delta_coord = 0;
                bool check_complex = true;
                for (size_t j = 0; j < l3.size(); j++) {
                    if (isnan(l3[j].real()) || isnan(l3[j].imag())) {
                        check_complex = false;
                        break;
                    }
                    delta_prev = std::max(delta_prev, std::abs(std::abs(l3[j]) - std::abs(l3_old[j])));
                    if (delta_prev > h || isnan(delta_prev)) {
                        check_complex = false;
                        break;
                    }

                    if (j >= 1)
                        delta_coord = std::max(delta_coord, std::abs(std::abs(l3[j]) - std::abs(l3[j-1])));
                    if (delta_coord > h || isnan(delta_coord)) {
                        check_complex = false;
                        break;
                    }
                }

                if ((l2[1] - l2[0]).norm() < h) {
                    std::cout << "PICKED SQUARES METHOD" << std::endl;
                    std::cout << "EIGENVEC" << lastVec[3] + (l2[1] * lastVec[2]) << std::endl;

                    return l2[1];

                } else if ((l1[2] - l1[1]).norm() < h) {
                    std::cout << "PICKED USUAL METHOD" << std::endl;
                    std::cout << "EIGENVEC" << lastVec[3] / lastVec[3].norm() << std::endl;

                    return l1[2];
                } else if (check_complex) {
                    std::cout << "PICKED COMPLEX METHOD" << std::endl;
                    std::cout << "EIGENVEC" << lastVec[3] + (l2[1] * lastVec[2]) << std::endl;

                    Vector result;
                    result(0) = l3[0].real();
                    return result;
                }

                u = u / u.norm();
                l3_old = l3;
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

std::vector<std::complex<double>> PowerIterations::findComplexL(Vector* lv) {
    //Vector l3(lv->length());
    std::vector<std::complex<double>> l3;
    for (int i = 0; i < lv->length(); i++) {
        double r_upper = lv[1](i)*lv[3](i) - lv[2](i)*lv[2](i);
        double r_lower = lv[0](i)*lv[2](i) - lv[1](i)*lv[1](i);
        double r = sqrt(r_upper / r_lower);
        double cos = (lv[3](i) + r*r * lv[1](i)) / (2 * r * lv[2](i));
        double sin = sqrt(1 - cos*cos);

        l3.emplace_back(r * cos, r * sin);
    }

    return l3;
}
