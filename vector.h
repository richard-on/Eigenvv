#ifndef MV_VECTOR_H
#define MV_VECTOR_H

#include <vector>
#include <complex>

class Vector {
public:
    explicit Vector(int len = 1);

    Vector(int len, double val);

    Vector(int len, double* data);

    Vector(std::initializer_list<double> list);

    Vector(const Vector& other);

    Vector(Vector&& other) noexcept;

    Vector& operator = (const Vector& other);

    Vector& operator = (Vector&& other) noexcept;

    double& operator () (int idx) const;

    bool operator == (const Vector& other) const;

    bool operator != (const Vector& other) const;

    Vector operator + (const Vector& a) const;

    Vector operator - (const Vector& a) const;

    Vector operator / (const Vector& other) const;

    Vector operator / (double other) const;

    Vector operator * (const Vector& other) const;

    friend std::ostream& operator << (std::ostream& ostream, const Vector& v);


    double* getData() const;

    int length() const;

    double norm();

    std::vector<std::complex<double>> toVector();


    virtual ~Vector();

private:
    double* data{};
    int len{};
};


#endif //MV_VECTOR_H
