#ifndef MV_02_NONLINEAR_H
#define MV_02_NONLINEAR_H

#include <vector>

struct Interval {
    double start;
    double end;
};

struct Solution {
    double root;
    Interval interval;
    int iterations;
};

class NonLinear {
public:
    NonLinear(double (*func)(double), double intervalStart, double intervalEnd, double step = 1e-1);

    Solution bisect(Interval &interval, double eps = 1e-4);

    Solution newton(Interval &interval, double eps = 1e-14, double maxIteration = 10000);

    std::vector<Solution> solve(double eps = 1e-14, double maxIteration = 1000);

    double df(double x, double eps = 1e-8);

    double d2f(double x, double eps = 1e-8);

    [[nodiscard]] const std::vector<Interval> &getRootIntervals() const;

private:
    double (*func)(double);
    Interval initialInterval{};
    std::vector<Interval> rootIntervals{};
};


#endif //MV_02_NONLINEAR_H
