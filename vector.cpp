#include <iostream>
#include <iomanip>
#include <climits>
#include <complex>

#include "vector.h"

Vector::Vector(int len) {
    if(len < 1) {
        throw std::invalid_argument("incorrect Vector length");
    }

    this->len = len;
    this->data = new double[len]{0};
}

Vector::Vector(int len, double val) {
    if(len < 1) {
        throw std::invalid_argument("incorrect Vector length");
    }

    this->len = len;
    this->data = new double[len];
    for(int i = 0; i < len; i++) {
        this->data[i] = val;
    }
}

Vector::Vector(int len, double* data) {
    if(len < 1) {
        throw std::invalid_argument("incorrect Vector length");
    }

    this->len = len;
    for(int i = 0; i < len; i++) {
        this->data[i] = data[i];
    }
}

Vector::Vector(std::initializer_list<double> list)
    :len{static_cast<int>(list.size())}, data{new double[static_cast<int>(list.size())]}{

    std::copy(list.begin(), list.end(), data);
}

Vector::Vector(const Vector &other) {
    this->len = other.len;
    this->data = new double[len]{0};
    for(int i = 0; i < len; i++) {
        this->data[i] = other.data[i];
    }
}

Vector::Vector(Vector &&other) noexcept
    : data(other.data), len(other.len)
    {

    other.data = nullptr;
    other.len = 0;
}


Vector& Vector::operator=(const Vector &other) {
    if(this != &other){
        delete[] this->data;

        this->len = other.len;
        this->data = new double[len];
        for(int i = 0; i < len; i++) {
            this->data[i] = other.data[i];
        }
    }

    return *this;
}

Vector& Vector::operator=(Vector &&other) noexcept {
    if(this != &other){
        delete[] data;

        this->len = other.len;
        this->data = new double[len];
        for(int i = 0; i < len; i++) {
            this->data[i] = other.data[i];
        }

        other.len = 0;
        other.data = nullptr;
    }

    return *this;
}

double &Vector::operator()(int i) const {
    return this->data[i];
}

bool Vector::operator==(const Vector &other) const {
    if(this->len == other.len){
        for(int i = 0; i < len; i++) {
            if(this->data[i] != other.data[i]) {
                return false;
            }
        }
        return true;
    }
    return false;
}

bool Vector::operator!=(const Vector &other) const {
    if(this == &other){
        return false;
    }
    return true;
}

Vector Vector::operator+(const Vector &other) const {
    if(this->len != other.len) {
        throw std::invalid_argument("error trying to compare vectors of different length");
    }

    Vector v(len);
    for (int i = 0; i < len; i++) {
        v.data[i] = this->data[i] + other.data[i];
    }

    return v;
}

Vector Vector::operator-(const Vector& other) const {
    if(this->len != other.len) {
        throw std::invalid_argument("error trying to compare vectors of different length");
    }

    Vector v(len);
    for (int i = 0; i < len; i++) {
        v.data[i] = this->data[i] - other.data[i];
    }

    return v;
}

Vector Vector::operator/(const Vector& other) const {
    if(this->len != other.len) {
        throw std::invalid_argument("error trying to compare vectors of different length");
    }

    Vector v(len);
    for (int i = 0; i < len; i++) {
        v.data[i] = this->data[i] / other.data[i];
    }

    return v;
}

Vector Vector::operator/(double other) const {
    Vector v(len);
    for (int i = 0; i < len; i++) {
        v.data[i] = this->data[i] / other;
    }

    return v;
}

Vector Vector::operator*(const Vector& other) const {
    Vector v(len);
    for (int i = 0; i < len; i++) {
        v(i) = this->data[i] * other.data[i];
    }

    return v;
}

std::ostream &operator<<(std::ostream &ostream, const Vector &vector) {
    ostream << std::endl;
    for (int i = 0; i < vector.length(); i++) {
        if(i == 0) {
            ostream << "/";
        }
        else if(i == vector.length() - 1) {
            ostream << "\\";
        }
        else{
            ostream << "|";
        }
        ostream << std::setw(8) << std::setprecision(3) << vector(i);
        if(i == 0) {
            ostream << "\\";
        }
        else if(i == vector.length() - 1) {
            ostream<< "/";
        }
        else{
            ostream << "|";
        }

        ostream << std::endl;
    }

    return ostream;
}


double *Vector::getData() const {
    return data;
}

int Vector::length() const {
    return len;
}

double Vector::norm() {
    double norm = 0;

    for (int i = 0; i < len; i++) {
        norm += std::pow(this->data[i], 2);
    }

    return std::sqrt(norm);
}


Vector::~Vector() {
    delete[] data;
}

std::vector<std::complex<double>> Vector::toVector() {
    std::vector<std::complex<double>> vec(len);

    for(int i = 0; i < len; i++) {
        vec[i] = data[i];
    }

    return vec;
}
