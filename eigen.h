#ifndef MV_02_EIGEN_H
#define MV_02_EIGEN_H

#include "vector.h"
#include "matrix.h"

enum class EigenKind {
    UsualMethod, SquareMethod, ComplexMethod, QRAlgorithm, NoValues
};

struct Eigenvv {
    std::complex<double> value;
    std::vector<std::complex<double>> vector;
};

class eigen {
public:
    eigen() = delete;

    static eigen powerIteration(const Matrix &matrix, Vector initialGuess = Vector(1),
                                double eps = 10e-8, int maxIteration = 20000);

    static eigen QR(Matrix m, double eps = 10e-8, int maxIteration = 20000);

    static std::string complexToString(std::complex<double> c);

    friend std::ostream& operator << (std::ostream& ostream, const eigen& eigen);

    static void normalize(std::vector<std::complex<double>> &cv);

    [[nodiscard]] EigenKind getEigenKind() const;

    [[nodiscard]] std::vector<std::complex<double>> getEigenValues() const;

    [[nodiscard]] std::vector<std::vector<std::complex<double>>> getEigenVectors() const;

    [[nodiscard]] int getNumberOfValues() const;

    [[nodiscard]] int getIterations() const;

private:
    EigenKind eigenKind;
    std::vector<Eigenvv> vv;
    int n;
    int iterations;

    eigen(EigenKind kind, std::vector<Eigenvv> vec, int n, int iter);
    static std::vector<std::complex<double>> findComplexLambda(std::vector<Vector> lv);
};


#endif //MV_02_EIGEN_H
