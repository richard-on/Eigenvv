#ifndef MV_02_EIGEN_H
#define MV_02_EIGEN_H

#include "vector.h"
#include "matrix.h"

enum class EigenKind {
    UsualMethod, SquareMethod, ComplexMethod, LeastSquares, NoVectors
};

struct Eigenvv {
    std::complex<double> value;
    std::vector<std::complex<double>> vector;
};

class Eigen {
public:
    Eigen() = delete;

    static Eigen powerIteration(const Matrix &matrix, Vector initialGuess, double h = 10e-8, int maxIteration = 10000);

    static Eigen QR(Matrix m, double h = 10e-8, int maxIteration = 10000);
    std::vector<std::complex<double>> ExtractEigenValuesQr(Matrix<double> m, double eps);
    static void RotateRow(Matrix &matrix, double c, double s, int master, int slave);
    static void RotateColumn(Matrix &matrix, double c, double s, int master, int slave);

    static std::string complexToString(std::complex<double> c);

    friend std::ostream& operator << (std::ostream& ostream, const Eigen& eigen);

private:
    EigenKind eigenKind;
    std::vector<Eigenvv> vv;
    int n;
    int iterations;

    Eigen(EigenKind kind, std::vector<Eigenvv> vector1, int i, int i1);
    static std::vector<std::complex<double>> findComplexLambda(Vector* lv);

};


#endif //MV_02_EIGEN_H
