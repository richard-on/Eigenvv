#ifndef MV_MATRIX_H
#define MV_MATRIX_H

#include <iosfwd>
#include <climits>
#include <vector>

#include "vector.h"

struct Sor{
    int step;
    std::vector<double> residual;
    Vector x;
};

class Matrix {
public:
    explicit Matrix(int len = 2);

    Matrix(int len, double** data);

    Matrix(std::initializer_list<std::initializer_list<double>> list);

    Matrix(const Matrix& other);

    Matrix(Matrix&& other) noexcept;

    Matrix& operator = (const Matrix& other);

    Matrix& operator = (Matrix &&other) noexcept;

    double& operator () (int row, int col) const;

    bool operator == (const Matrix &other) const;

    bool operator != (const Matrix &other) const;

    Matrix operator + (const Matrix &other) const;

    Matrix operator - (const Matrix &other) const;

    //Multiply matrix by another matrix, result is matrix
    Matrix operator * (const Matrix &other) const;

    //Multiply matrix by vector, result is vector
    Vector operator * (const Vector &other) const;

    friend std::ostream& operator << (std::ostream& ostream, const Matrix& matrix);

    static Matrix identity(int len);

    double** getData() const;

    int length() const;

    double norm();

    double condition();

    Matrix inverse();

    Matrix transpose() const;

    Vector solveGauss(const Vector& b);

    Vector solveGaussCol(const Vector& b);

    Sor solveSOR(const Vector& b, double param, double eps = 1e-8);


    virtual ~Matrix();

private:
    double** data{};
    int len{};
};


#endif //MV_MATRIX_H
