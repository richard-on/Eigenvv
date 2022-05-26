#include <complex>
#include <iostream>
#include <utility>

#include "eigen.h"

eigen::eigen(EigenKind eigenKind, std::vector<Eigenvv> vec, int n, int iterations) {
    this->eigenKind = eigenKind;
    this->vv = std::move(vec);
    this->n = n;
    this->iterations = iterations;
}

eigen eigen::powerIteration(const Matrix &matrix, Vector u0, double eps, int maxIteration) {
    if (u0.length() != matrix.length()) {
        u0 = Vector(matrix.length(), 1);
    }

    Vector v(u0.length());
    Vector u(u0.length());

    std::vector<Vector> vec(4);
    std::vector<Vector> l1(3);
    std::vector<Vector> l2(2);

    std::vector<std::complex<double>> curL3(u0.length(), 0);
    std::vector<std::complex<double>> prevL3(u0.length(), 0);

    for (int i = 0; i < maxIteration; i++) {
        u = (matrix * u0);
        switch (i%4) {
            case 0:
                vec[0] = u;
                break;
            case 1:
                vec[1] = u;
                break;
            case 2:
                vec[2] = u;
                break;
            case 3:
                vec[3] = u;

                l1[0] = vec[1] / vec[0];
                l1[1] = vec[2] / vec[1];
                l1[2] = vec[3] / vec[2];

                l2[0] = vec[2] / vec[0];
                l2[1] = vec[3] / vec[1];

                for(int j = 0; j < l2[0].length(); j++) {
                    l2[0](j) = std::sqrt(std::abs(l2[0](j)));
                    l2[1](j) = std::sqrt(std::abs(l2[1](j)));
                }

                curL3 = findComplexLambda(vec);

                double prevDelta = 0;
                double curDelta = 0;
                bool check_complex = true;
                for (int j = 0; j < curL3.size(); j++) {
                    if (std::isnan(curL3[j].real()) || std::isnan(curL3[j].imag())) {
                        check_complex = false;
                        break;
                    }

                    prevDelta = std::max(prevDelta, std::abs(std::abs(curL3[j]) - std::abs(prevL3[j])));
                    if (prevDelta > eps || std::isnan(prevDelta)) {
                        check_complex = false;
                        break;
                    }

                    if (j >= 1) {
                        curDelta = std::max(curDelta, std::abs(std::abs(curL3[j]) - std::abs(curL3[j - 1])));
                    }

                    if (curDelta > eps || std::isnan(curDelta)) {
                        check_complex = false;
                        break;
                    }
                }

                if ((l1[2] - l1[1]).norm() < eps) {
                    std::vector<Eigenvv> eigenVec;
                    eigenVec.push_back(Eigenvv{l1[2](0), vec[3].toVector()});

                    eigen::normalize(eigenVec[0].vector);

                    return {EigenKind::UsualMethod,
                            eigenVec,
                            static_cast<int>(eigenVec.size()),
                            i
                    };

                } else if ((l2[1] - l2[0]).norm() < eps) {
                    std::vector<Eigenvv> eigenVec;
                    eigenVec.push_back(Eigenvv{l2[1](0), vec[3].toVector()});
                    eigenVec.push_back(Eigenvv{-l2[1](0), vec[3].toVector()});

                    for (int j = 0; j < vec[3].length(); j++) {
                        eigenVec[0].vector[j] = vec[3].toVector()[j] + eigenVec[0].value * vec[1].toVector()[j];
                        eigenVec[1].vector[j] = vec[3].toVector()[j] + eigenVec[1].value * vec[1].toVector()[j];
                    }

                    eigen::normalize(eigenVec[0].vector);
                    eigen::normalize(eigenVec[1].vector);
                    return {EigenKind::SquareMethod,
                            eigenVec,
                            static_cast<int>(eigenVec.size()),
                            i
                    };

                } else if (check_complex) {
                    std::vector<Eigenvv> eigenVec;
                    eigenVec.push_back(Eigenvv{curL3[0], vec[0].toVector()});
                    eigenVec.push_back(Eigenvv{{curL3[0].real(), -curL3[0].imag()}, vec[3].toVector()});

                    for (int j = 0; j < vec[3].length(); j++) {
                        eigenVec[0].vector[j] = {vec[3].toVector()[j].real() - vec[2].toVector()[j].real() * eigenVec[0].value.real(),
                                                 vec[2].toVector()[j].real() * eigenVec[0].value.imag()};

                        eigenVec[1].vector[j] = {vec[3].toVector()[j].real() - vec[2].toVector()[j].real() * eigenVec[1].value.real(),
                                                 vec[2].toVector()[j].real() * eigenVec[1].value.imag()};

                    }
                    eigen::normalize(eigenVec[0].vector);
                    eigen::normalize(eigenVec[1].vector);

                    return {EigenKind::ComplexMethod,
                            eigenVec,
                            static_cast<int>(eigenVec.size()),
                            i
                    };
                }

                u = u / u.norm();
                prevL3 = curL3;
                break;
        }
        u0 = u;

    }

    return {EigenKind::NoValues,std::vector<Eigenvv>{},0,maxIteration};
}

std::vector<std::complex<double>> eigen::findComplexLambda(std::vector<Vector> lv) {
    std::vector<std::complex<double>> l3;
    for (int i = 0; i < lv.size(); i++) {
        double r_upper = lv[1](i)*lv[3](i) - lv[2](i)*lv[2](i);
        double r_lower = lv[0](i)*lv[2](i) - lv[1](i)*lv[1](i);
        double r = sqrt(r_upper / r_lower);
        double cos = (lv[3](i) + r*r * lv[1](i)) / (2 * r * lv[2](i));
        double sin = sqrt(1 - cos*cos);

        l3.emplace_back(r * cos, r * sin);
    }

    return l3;
}

void eigen::normalize(std::vector<std::complex<double>> &cv) {
    double norm = 0;
    for (auto & i : cv) {
        norm += std::abs(i * i);
    }

    for (auto & i : cv) {
        i /= sqrt(norm);
    }

}

eigen eigen::QR(Matrix m, double eps, int maxIteration) {
    m = m.hessenberg();

    std::vector<std::complex<double>> curEigenv;
    int prevSize = 0;
    int requiredZeroCount = (m.length() - 1) / 2;
    bool hasEigenv = false;

    int iter = 0;
    while (iter < maxIteration) {
        Matrix q = Matrix::identity(m.length());

        for (int i = 0; i < m.length() - 1; i++) {
            double csRoot = std::sqrt(std::pow(m(i, i), 2) + std::pow(m(i + 1, i), 2));
            double cos = m(i, i) / csRoot;
            double sin = m(i + 1, i) / csRoot;

            for (int j = 0; j < m.length(); j++) {
                double tempM = m(i, j);
                double tempQ = q(j, i);

                m(i, j) = cos * m(i, j) + sin * m(i + 1, j);
                m(i + 1, j) = -sin * tempM + cos * m(i + 1, j);

                q(j, i) = cos * q(j, i) + sin * q(j, i + 1);
                q(j, i + 1) = -sin * tempQ + cos * q(j, i + 1);
            }
        }
        m = m * q;

        if (hasEigenv) {
            std::vector<int> value_type = std::vector<int>(m.length() + 1, 0);
            curEigenv = {};

            for (int i = 0; i < m.length() - 1; i++) {
                if (std::abs(m(i+1, i)) > eps) {
                    value_type[i]++;
                    value_type[i + 1]++;
                    if (value_type[i] > 1) {
                        curEigenv = {};
                    }
                }
            }

            for (int i = 0; i < m.length(); i++) {
                if (value_type[i] == 0) {
                    curEigenv.emplace_back(m(i, i));

                } else if (value_type[i] == 1 && value_type[i+1] == 1) {
                    double b = m(i, i) + m(i + 1, i + 1);
                    double d = std::pow(b, 2) - 4 * (m(i, i)
                            * m(i + 1, i + 1) - m(i + 1, i) * m(i, i + 1));

                    if (d >= 0) {
                        curEigenv = {};
                    } else {
                        d = sqrt(-d);
                    }

                    double real = b / 2;
                    double complex = d / 2;
                    curEigenv.emplace_back(real, complex);
                    curEigenv.emplace_back(real, -complex);
                    i++;
                }
            }

            if (curEigenv.size() == m.length() && prevSize == m.length()) {
                break;
            }
            prevSize = curEigenv.size();

        } else {
            int subDiagZeroCount = 0;
            for (int i = 0; i < m.length() - 1; i++) {
                if (std::abs(q(i + 1, i)) < eps) {
                    subDiagZeroCount++;
                }
            }
            if (subDiagZeroCount >= requiredZeroCount) {
                hasEigenv = true;
            }

        }
        iter++;

    }

    std::vector<Eigenvv> vv(curEigenv.size());
    for (int i = 0; i < curEigenv.size(); i++) {
        vv[i].value = curEigenv[i];
        vv[i].vector = {};
    }
    return {
        EigenKind::QRAlgorithm,
        vv,
        static_cast<int>(curEigenv.size()),
        iter
    };
}

std::ostream &operator<<(std::ostream &ostream, const eigen &eigen) {

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
        case EigenKind::QRAlgorithm:
            ostream << "QR algorithm" << std::endl;
            break;
        case EigenKind::NoValues:
            ostream << "No values" << std::endl;
            return ostream;
        default:
            ostream << "No values" << std::endl;
            return ostream;
    }

    ostream << "Number of EigenValues found: " << eigen.n << std::endl;
    ostream << "Number of iterations: " << eigen.iterations << std::endl;

    for (int i = 0; i < eigen.vv.size(); i++) {
        ostream << "Eigen Value " << i + 1 << ": ";
        ostream << eigen::complexToString(eigen.vv[i].value) << std::endl;

        if (!eigen.vv[i].vector.empty()) {
            ostream << "Corresponding Eigen Vector: (";
        }

        int k = 0;
        for (auto j : eigen.vv[i].vector) {
            ostream << eigen::complexToString(j);
            if (k < eigen.vv[i].vector.size() - 1) {
                ostream << ", ";
            }
            k++;
        }
        if (!eigen.vv[i].vector.empty()) {
            ostream << ");" << std::endl;
        }

    }

    return ostream;
}

std::string eigen::complexToString(std::complex<double> c) {

    if (c.imag() == 0) {
        return std::to_string(c.real());
    } else if (c.imag() < 0) {
       return std::to_string(c.real()) + "-i*" + std::to_string(std::abs(c.imag()));
    } else {
        return std::to_string(c.real()) + "+i*" + std::to_string(c.imag());
    }
}

EigenKind eigen::getEigenKind() const {
    return eigenKind;
}

std::vector<std::complex<double>> eigen::getEigenValues() const {
    std::vector<std::complex<double>> values(vv.size());
    for (int i = 0; i < vv.size(); i++) {
        values[i] = vv[i].value;
    }

    return values;
}

std::vector<std::vector<std::complex<double>>> eigen::getEigenVectors() const {
    std::vector<std::vector<std::complex<double>>> vectors(vv.size());
    for (int i = 0; i < vv.size(); i++) {
        vectors[i] = vv[i].vector;
    }

    return vectors;
}

int eigen::getNumberOfValues() const {
    return n;
}

int eigen::getIterations() const {
    return iterations;
}
