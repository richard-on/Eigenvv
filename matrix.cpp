#include <iostream>
#include <iomanip>
#include <climits>

#include "matrix.h"

Matrix::Matrix(int len) {
    if(len < 2) {
        throw std::invalid_argument("incorrect Matrix dimensions");
    }

    this->len = len;
    this->data = new double*[len];
    for(int i = 0; i < len; i++) {
        this->data[i] = new double[len]{0};
    }
}

Matrix::Matrix(int len, double** data) {
    if(len < 2) {
        throw std::invalid_argument("incorrect Matrix dimensions");
    }

    this->len = len;
    this->data = new double*[len];
    for(int i = 0; i < len; i++) {
        this->data[i] = new double[len];
        for(int j = 0; j < len; j++) {
            this->data[i][j] = data[i][j];
        }
    }
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list)
    : len{static_cast<int>(list.size())} {

    for (auto val : list) {
        if (val.size() != len) {
            throw std::invalid_argument("not a square matrix");
        }
    }

    data = new double*[len];
    int i = 0;
    for (auto val : list) {
        data[i] = new double[len];
        std::copy(val.begin(), val.end(), data[i]);
        i++;
    }
}

Matrix::Matrix(const Matrix &other) {
    this->len = other.len;

    this->data = new double*[len];
    for(int i = 0; i < len; i++) {
        this->data[i] = new double[len];
        for(int j = 0; j < len; j++) {
            this->data[i][j] = other(i, j);
        }
    }
}

Matrix::Matrix(Matrix &&other) noexcept
    : data(other.data), len(other.len)
    {

    other.data = nullptr;
    other.len = 0;
}


Matrix &Matrix::operator=(const Matrix& other) {
    if(this != &other){
        for(int i = 0; i < len; i++) {
            delete[] data[i];
        }
        delete[] data;

        this->len = other.len;
        this->data = new double*[len];
        for(int i = 0; i < len; i++) {
            this->data[i] = new double[len];
            for(int j = 0; j < len; j++) {
                this->data[i][j] = other.data[i][j];
            }
        }
    }

    return *this;
}

Matrix &Matrix::operator=(Matrix &&other) noexcept {
    if (this != &other) {
        for (int i = 0; i < len; i++) {
            delete[] data[i];
        }
        delete[] data;

        this->len = other.len;
        this->data = new double*[len];
        for (int i = 0; i < len; i++) {
            this->data[i] = new double[len];
            for (int j = 0; j < len; j++) {
                this->data[i][j] = other.data[i][j];
            }
        }

        other.data = nullptr;
        other.len = 0;
    }

    return *this;
}

double& Matrix::operator()(int row, int col) const {
    return this->data[row][col];
}

bool Matrix::operator==(const Matrix &other) const {
    if(this->len == other.len){
        for(int i = 0; i < len; i++) {
            for(int j = 0; j < len; j++) {
                if(this->data[i][j] != other.data[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }
    return false;
}

bool Matrix::operator!=(const Matrix &other) const {
    if(this == &other){
        return false;
    }
    return true;
}

Matrix Matrix::operator+(const Matrix &other) const {
    if(this->len != other.len) {
        throw std::invalid_argument("error trying to compare Matrix of different length");
    }

    Matrix m = Matrix(len);
    for(int i = 0; i < len; i++) {
        for(int j = 0; j < len; j++) {
            m.data[i][j] = this->data[i][j] + other.data[i][j];
        }
    }

    return m;
}

Matrix Matrix::operator-(const Matrix &other) const {
    if(this->len != other.len) {
        throw std::invalid_argument("error trying to compare Matrix of different length");
    }

    Matrix m = Matrix(len);
    for(int i = 0; i < len; i++) {
        for(int j = 0; j < len; j++) {
            m.data[i][j] = this->data[i][j] - other.data[i][j];
        }
    }

    return m;
}

Matrix Matrix::operator*(const Matrix &other) const {
    if(this->len != other.len) {
        throw std::invalid_argument("error trying to compare Matrix of different length");
    }

    Matrix m = Matrix(len);
    for(int i = 0; i < len; i++) {
        for(int j = 0; j < len; j++) {
            for(int k = 0; k < len; k++) {
                m.data[i][j] += this->data[i][k] * other.data[k][j];
            }
        }
    }

    return m;
}

Vector Matrix::operator*(const Vector &other) const {
    if(other.length() != this->len) {
        throw std::invalid_argument("error trying to compare Matrix of different length");
    }

    Vector v(len);
    for(int i = 0; i < v.length(); i++) {
        for(int k = 0; k < v.length(); k++) {
            v(i) += data[i][k] * other(k);
        }
    }

    return v;
}

std::ostream &operator<<(std::ostream &ostream, const Matrix &matrix) {
    ostream << std::endl;
    for (int i = 0; i < matrix.length(); i++) {
        if(i == 0) {
            ostream << "/";
        }
        else if(i == matrix.length() - 1) {
            ostream << "\\";
        }
        else{
            ostream << "|";
        }

        for (int j = 0; j < matrix.length(); j++) {

            ostream << std::setw(10) << std::setprecision(3) << matrix(i, j) << " ";
        }
        if(i == 0) {
            ostream << "\\";
        }
        else if(i == matrix.length() - 1) {
            ostream << "/";
        }
        else{
            ostream << "|";
        }
        ostream << std::endl;

    }

    return ostream;
}


double **Matrix::getData() const {
    return data;
}

int Matrix::length() const {
    return len;
}

Matrix Matrix::identity(int len) {
    Matrix m(len);
    for(int i = 0; i < len; i++) {
        m.data[i][i] = 1;
    }

    return m;
}

Matrix Matrix::hessenberg() {
    Matrix m = *this;

    for(int i = len - 2; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {

            if (m.data[i + 1][i] == 0) {
                for (int swap = 0; swap < len; swap++) {
                    std::swap(m.data[swap][i], m.data[swap][j]);
                    std::swap(m.data[i][swap], m.data[j][swap]);
                }
            }

            double mult = -(m.data[i+1][j] / m.data[i+1][i]);
            for (int swap = 0; swap < len; swap++) {
                m.data[swap][j] += m.data[swap][i] * mult;
                m.data[i][swap] -= m.data[j][swap] * mult;
            }

        }
    }

    return m;
}

double Matrix::norm() {
    double norm = INT_MIN;
    double columnSum = 0;

    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            columnSum += std::abs(this->data[i][j]);
        }
        if(columnSum > norm) {
            norm = columnSum;
        }
        columnSum = 0;
    }

    return norm;
}

double Matrix::condition() {
    double norm = this->norm();
    double invNorm = this->inverse().norm();

    return norm * invNorm;
}

Matrix Matrix::inverse() {
    auto res = Matrix(len);
    auto** m = new double*[len];
    for (int i = 0; i < len; i++) {
        m[i] = new double[len * 2]{0};
        for (int j = 0; j < 2 * len; j++) {
            if(j < len) {
                m[i][j] = this->data[i][j];
            }
            if (j == (i + len))
                m[i][j] = 1;
        }
    }

    for (int i = len-1; i > 0; i--) {
        if (m[i - 1][0] < m[i][0]) {
            double* temp = m[i];
            m[i] = m[i - 1];
            m[i - 1] = temp;
        }
    }

    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            if (j != i) {
                double d = m[j][i] / m[i][i];
                for (int k = 0; k < 2 * len; k++) {
                    m[j][k] -= m[i][k] * d;
                }
            }
        }
    }

    for (int i = 0; i < len; i++) {
        double d = m[i][i];
        for (int j = 0; j < 2 * len; j++) {
            m[i][j] = m[i][j] / d;
        }
    }

    for (int i = 0; i < len; i++) {
        int k = 0;
        for (int j = len; j < 2 * len; j++) {
            res(i, k) = m[i][j];
            k++;
        }
    }

    for (int i = 0; i < len; i++) {
        delete[] m[i];
    }
    delete[] m;

    return res;
}

Matrix Matrix::transpose() const {
    Matrix t = Matrix(this->len);
    for (int i = 0; i < this->len; i++) {
        for (int j = 0; j < this->len; j++) {
            t(i, j) = (*this)(j, i);
        }
    }

    return t;
}

Vector Matrix::solveGauss(const Vector &b) {
    Vector x(b.length());

    auto** res = new double*[this->len];
    for (int i = 0; i < this->len; i++) {
        res[i] = new double[this->len + 1]{0};
        for(int j = 0; j < this->len + 1; j++) {
            if(j == this->len) {
                res[i][j] = b(i);
            }
            else {
                res[i][j] = this->data[i][j];
            }
        }
    }

    for (int i = 0; i < this->len - 1; i++) {
        for (int j = i + 1; j < this->len; j++) {
            double d = res[j][i]/res[i][i];
            for(int k = 0; k < this->len + 1; k++) {
                res[j][k] = res[j][k]-d*res[i][k];
            }
        }
    }

    for(int i = this->len - 1; i >= 0; i--) {
        x(i) = res[i][this->len];
        for(int j = i + 1; j < this->len; j++) {
            if(i != j) {
                x(i) = x(i) - (res[i][j] * x(j));
            }
        }
        x(i) = x(i) / res[i][i];
    }

    for(int i = 0; i < len; i++){
        delete[] res[i];
    }
    delete[] res;

    return x;
}

Vector Matrix::solveGaussCol(const Vector &b) {
    Vector x(b.length());

    double aa, bb;

    auto** res = new double*[this->len];
    for (int i = 0; i < this->len; i++) {
        res[i] = new double[this->len + 1]{0};
        for(int j = 0; j < this->len + 1; j++) {
            if(j == this->len) {
                res[i][j] = b(i);
            }
            else {
                res[i][j] = this->data[i][j];
            }
        }
    }

    for (int i = 0, k = 0; k < this->len; k++) {
        aa = std::abs(res[k][k]);
        i = k;
        for(int m = k + 1; m < this->len; m++) {
            if (std::abs(res[m][k]) > aa) {
                i = m;
                aa = std::abs(res[m][k]);
            }
        }

        if (i != k) {
            for (int j = k; j < this->len + 1; j++) {
                bb = res[k][j];
                res[k][j] = res[i][j];
                res[i][j] = bb;
            }
        }

        aa = res[k][k];
        res[k][k] = 1;
        for (int j = k + 1; j < this->len + 1; j++) {
            res[k][j] = res[k][j] / aa;
        }
        for (i = k + 1; i < this->len; i++) {
            bb = res[i][k];
            res[i][k] = 0;
            if (bb != 0) {
                for (int j = k + 1; j < this->len + 1; j++) {
                    res[i][j] = res[i][j] - bb * res[k][j];
                }
            }
        }
    }

    for(int i = this->len - 1; i >= 0; i--) {
        x(i) = res[i][this->len];
        for(int j = i + 1; j < this->len; j++) {
            if(i != j) {
                x(i) = x(i) - (res[i][j] * x(j));
            }
        }
        x(i) = x(i) / res[i][i];
    }

    for(int i = 0; i < len; i++){
        delete[] res[i];
    }
    delete[] res;

    return x;
}

Sor Matrix::solveSOR(const Vector &b, double param, double eps) {
    Vector x(b.length());
    std::vector<double> residual;

    int step = 0;
    double norm;
    double res = (*this * x - b).norm();

    while(res > eps) {
        for(int i = 0; i < this->len; i++) {
            norm = 0;
            for(int j = 0; j < this->len; j++) {
                if(j != i) {
                    norm += this->data[i][j]*x(j);
                }
            }
            x(i) = (1 - param) * x(i) + (param / this->data[i][i]) * (b(i) - norm);
        }
        res = (*this * x - b).norm();
        residual.push_back(res);
        step++;
        if (step > 50) {
            break;
        }
    }

    Sor sorSolution{step, residual, x};

    return sorSolution;
}


Matrix::~Matrix() {
    for (int i = 0; i < len; i++) {
        delete[] data[i];
    }
    delete[] data;
}

