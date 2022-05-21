#ifndef MV_02_NONLINEAR_H
#define MV_02_NONLINEAR_H

struct Interval {
    double start;
    double end;
};

class NonLinear {
public:
    NonLinear(double (*func)(double), double intervalStart, double intervalEnd, double step = 10e-2);

    Interval bisect(int rootPos, double h = 10e-3);

    double newton(int rootPos, double h = 10e-8, double maxIteration = 50000);

    std::vector<double> solve(double h = 10e-8, double maxIteration = 50000);

    double derivative(double x);

    const Interval &getInterval() const;

    int getRootsNum() const;

    const std::vector<Interval> &getRoots() const;

private:
    double (*func)(double);
    Interval interval{};
    int rootsNum;
    std::vector<Interval> roots;
};


#endif //MV_02_NONLINEAR_H
