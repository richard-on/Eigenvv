#ifndef MV_02_POWERITERATIONS_H
#define MV_02_POWERITERATIONS_H


#include "vector.h"
#include "matrix.h"

class PowerIterations {
public:
    static Vector powerIteration(const Matrix &matrix, Vector initialGuess, double h, int maxIteration = 10000);

private:
    static std::vector<std::complex<double>> findComplexL(Vector* lv);

};


#endif //MV_02_POWERITERATIONS_H
