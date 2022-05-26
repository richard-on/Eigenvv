#include <vector>
#include <cmath>

#include "nonLinear.h"

NonLinear::NonLinear(double (*func)(double), double intervalStart, double intervalEnd, double step) {
    std::vector<Interval> intervals;
    this->func = func;
    this->initialInterval.start = intervalStart;
    this->initialInterval.end = intervalEnd;
    int rootNum = 0;

    double x1 = initialInterval.start;
    double x2 = x1 + step;
    double y1 = func(x1);

    while (x2 < initialInterval.end) {
        double y2 = func(x2);
        if (y1 <= 0 && y2 >= 0 || y1 >= 0 && y2 <= 0) {
            intervals.push_back(Interval{x1, x2});
            rootNum++;
        }

        x1 = x2;
        x2 = x1 + step;
        y1 = y2;
    }

    this->rootIntervals = intervals;
}

Solution NonLinear::bisect(Interval &interval, double eps) {
    if (func(interval.start) == 0) {
        return Solution{interval.start, interval, 0};
    }
    if (func(interval.end) == 0) {
        return Solution{interval.end, interval, 0};
    }

    double x0 = interval.start;
    double x = interval.end;
    double xi;
    int iter = 0;

    while (std::abs(x0 - x) > eps) {
        xi = x + (x0 - x) / 2;
        if(func(x) < 0 && func(xi) > 0 || func(x) > 0 && func(xi) < 0) {
            x0 = xi;
        } else {
            x = xi;
        }

        iter++;
    }

    interval.start = x0;
    interval.end = x;

    return Solution{(x - x0)/2, interval, iter};
}

Solution NonLinear::newton(Interval &interval, double eps, double maxIteration) {
    double x0 = interval.start;
    int iter = 0;
    if (func(x0) >= 0 && d2f(x0) >= 0 || func(x0) < 0 && d2f(x0) < 0) {
        x0 = interval.start;
    } else {
        x0 = interval.end;
    }

    double x = x0 - func(x0) / df(x0);
    while (std::abs(x0 - x) > eps) {
        x0 = x;
        x = x - func(x0) / df(x0);
        iter++;
        if (iter > maxIteration) {
            break;
        }
    }

    return Solution{x, {x, x}, iter};
}

std::vector<Solution> NonLinear::solve(double eps, double maxIteration) {
    std::vector<Solution> solutions;

    for (auto& rootInterval : rootIntervals) {
        bisect(rootInterval, 10e-4);
        solutions.emplace_back(newton(rootInterval, eps, maxIteration));
    }

    return solutions;
}

double NonLinear::df(double x, double eps) {
    return (func(x + eps) - func(x - eps)) / (2 * eps);
}

double NonLinear::d2f(double x, double eps) {
    return (func(x + eps) - 2 * func(x) + func(x - eps)) / (eps * eps);
}

const std::vector<Interval> &NonLinear::getRootIntervals() const {
    return rootIntervals;
}
