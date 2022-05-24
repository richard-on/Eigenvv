#include <complex>
#include <iostream>
#include <array>

#include "eigen.h"

Eigen::Eigen(EigenKind eigenKind, std::vector<Eigenvv> vec, int n, int iterations) {
    this->eigenKind = eigenKind;
    this->vv = vec;
    this->n = n;
    this->iterations = iterations;
}

Eigen Eigen::powerIteration(const Matrix &a, Vector u0, double h, int maxIteration) {
    Vector v(u0.length());
    Vector u(u0.length());
    Vector um1(u.length());

    auto* lastVec = new Vector[4];
    auto* l1 = new Vector[3];
    auto* l2 = new Vector[2];
    std::vector<std::complex<double>> l3(u0.length(), 0);
    std::vector<std::complex<double>> l3_old(u0.length(), 0);
    for (int i = 0; i < maxIteration; i++) {
        u = (a * u0);
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
                for(int j = 0; j < l2[0].length(); j++) {
                    l2[0](j) = std::sqrt(std::abs(l2[0](j)));
                    l2[1](j) = std::sqrt(std::abs(l2[1](j)));
                }

                l3 = findComplexLambda(lastVec);

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

                    std::vector<Eigenvv> eigenVec;
                    eigenVec.emplace_back(l2[1](0), l2[1].toVector());
                    eigenVec.emplace_back(-l2[1](0), l2[1].toVector());

                    return {EigenKind::SquareMethod,
                            eigenVec,
                            static_cast<int>(eigenVec.size()),
                            i
                    };

                } else if ((l1[2] - l1[1]).norm() < h) {
                    std::cout << "PICKED USUAL METHOD" << std::endl;
                    std::cout << "EIGENVEC" << lastVec[3] / lastVec[3].norm() << std::endl;

                    std::vector<Eigenvv> eigenVec;
                    eigenVec.emplace_back(l1[2](0), l2[1].toVector());

                    return {EigenKind::UsualMethod,
                            eigenVec,
                            static_cast<int>(eigenVec.size()),
                            i
                    };

                } else if (check_complex) {
                    std::cout << "PICKED COMPLEX METHOD" << std::endl;
                    std::cout << "EIGENVEC" << lastVec[3] + (l2[1] * lastVec[2]) << std::endl;

                    std::vector<Eigenvv> eigenVec;
                    eigenVec.emplace_back(l3[0], l3);

                    return {EigenKind::ComplexMethod,
                            eigenVec,
                            static_cast<int>(eigenVec.size()),
                            i
                    };
                }

                u = u / u.norm();
                l3_old = l3;
                break;
        }

        u0 = u;

    }

    return {EigenKind::NoVectors,std::vector<Eigenvv>{},0,maxIteration};
}

std::vector<std::complex<double>> Eigen::findComplexLambda(Vector* lv) {
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

static std::vector<std::complex<double>> ExtractEigenValuesQr(Matrix<double> m, double eps) {
    enum EigenValueType : int {
        REAL = 0,
        COMPLEX = 1,
        NOT_READY = 2
    };
    std::vector<int> value_type = std::vector<int>(m.size_, 0);

    for (size_t i = 0; i < m.size_ - 1; i++) {
        if (std::abs(m.data[i+1][i]) > eps) {
            ++value_type[i];
            ++value_type[i + 1];
            if (value_type[i] == NOT_READY) {
                return {}; // Need more iterations
            }
        }
    }

    std::vector<std::complex<double>> eigen_values;
    for (size_t i = 0; i < m.size_; i++) {
        if (value_type[i] == REAL) {
            eigen_values.push_back(m.data[i][i]);
        } else if (value_type[i] == COMPLEX && value_type[i+1] == COMPLEX) {
            double b = m.data[i][i] + m.data[i + 1][i + 1];
            double d = b * b - 4 * (m.data[i][i] * m.data[i + 1][i + 1] - m.data[i + 1][i] * m.data[i][i + 1]);
            if (d >= 0) {
                return {}; // Need more iterations
            }
            double sqrt_d = sqrt(-d);

            double real = b / 2;
            double complex = sqrt_d / 2;
            eigen_values.emplace_back(real, complex);
            eigen_values.emplace_back(real, -complex);
            ++i; // skip one step, complex pair was found
        }
    }
    return eigen_values;
}

Eigen Eigen::QR(Matrix m, double h, int maxIteration) {
    size_t min_zeroes = (m.length() - 1) / 2;
    bool is_get_min_zeroes = false;
    std::vector<std::complex<double>> eig_vals_cur{};
    std::vector<std::complex<double>> eig_vals_prev{};
    Matrix q = Matrix::identity(m.length());


    for (int i = 0; i < m.length() - 1; i++) {
        double root = std::sqrt(std::pow(m(i, i), 2) + std::pow(m(i + 1, i), 2));
        double c = m(i, i) / root;
        double s = m(i + 1, i) / root;

        RotateRow(m, c, s, i, i + 1);
        RotateColumn(q, c, -s, i, i + 1);
    }
    m = m * q;

    if (is_get_min_zeroes) {
        eig_vals_cur = ExtractEigenValuesQr(m, h);
        if (eig_vals_cur.size() == m.length() && eig_vals_prev.size() == m.length()) {
            double delta = 0;
            for (size_t i = 0; i < eig_vals_cur.size(); i++) {
                std::max(delta, std::abs(std::abs(eig_vals_cur[i]) - std::abs(eig_vals_prev[i])));
            }
            if (delta < h) {
                return eig_vals_cur;
            }
        }
        eig_vals_prev = std::move(eig_vals_cur);
    } else {
        size_t zeroes_count = 0;
        for (size_t i = 0; i < m.length() - 1; i++) {
            if (std::abs(q(i + 1, i)) < h) {
                ++zeroes_count;
            }
        }
        is_get_min_zeroes = zeroes_count >= min_zeroes;
    }


}

void Eigen::RotateRow(Matrix &matrix, double c, double s, int master, int slave) {
    for (int i = 0; i < matrix.length(); i++) {
        double temp = matrix(master, i);
        matrix(master, i) = c * matrix(master, i) + s * matrix(slave, i);
        matrix(slave, i) = -s * temp + c * matrix(slave, i);
    }
}

void Eigen::RotateColumn(Matrix &matrix, double c, double s, int master, int slave) {
    for (int i = 0; i < matrix.length(); i++) {
        double temp = matrix(i, master);
        matrix(i, master) = c * matrix(i, master) + -s * matrix(i, slave);
        matrix(i, slave) = s * temp + c * matrix(i, slave);
    }
}


std::ostream &operator<<(std::ostream &ostream, const Eigen &eigen) {

    switch (eigen.eigenKind) {
        case EigenKind::UsualMethod:
            ostream << "Usual Method" << std::endl;
            break;
        case EigenKind::SquareMethod:
            ostream << "Square Method" << std::endl;
            break;
        case EigenKind::ComplexMethod:
            ostream << "Complex Method" << std::endl;
            break;
        case EigenKind::LeastSquares:
            ostream << "Least Squares Method" << std::endl;
            break;
        case EigenKind::NoVectors:
            ostream << "No vectors" << std::endl;
            return ostream;
        default:
            ostream << "No vectors" << std::endl;
            return ostream;
    }

    ostream << "Number of EigenValues found: " << eigen.n << std::endl;
    ostream << "Number of iterations: " << eigen.iterations << std::endl;

    for (int i = 0; i < eigen.vv.size(); i++) {
        ostream << "Eigen Value " << i + 1 << ": ";
        ostream << Eigen::complexToString(eigen.vv[i].value) << std::endl;
        ostream << "Corresponding Eigen Vector: (";

        int k = 0;
        for (auto j : eigen.vv[i].vector) {
            ostream << Eigen::complexToString(j);
            if (k < eigen.vv[i].vector.size() - 1) {
                ostream << ", ";
            }
            k++;
        }
        ostream << ");" << std::endl;
        ostream << std::endl;
    }

    return ostream;
}

std::string Eigen::complexToString(std::complex<double> c) {

    if (c.imag() == 0) {
        return std::to_string(c.real());
    } else if (c.imag() < 0) {
       return std::to_string(c.real()) + "-i*" + std::to_string(c.imag());
    } else {
        return std::to_string(c.real()) + "+i*" + std::to_string(c.imag());
    }
}
