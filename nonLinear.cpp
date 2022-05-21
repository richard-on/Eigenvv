#include <vector>
#include "nonLinear.h"

NonLinear::NonLinear(double (*func)(double), double intervalStart, double intervalEnd, double step) {
    std::vector<Interval> intervals;
    this->func = func;
    this->interval.start = intervalStart;
    this->interval.end = intervalEnd;
    int rootNum = 0;

    double x1 = interval.start;
    double x2 = x1 + step;
    double y1 = func(x1);

    while (x2 < interval.end) {
        double y2 = func(x2);
        if (y1 * y2 < 0) {
            intervals.emplace_back(Interval{x1, x2});
            rootNum++;
        }

        x1 = x2;
        x2 = x1 + step;
        y1 = y2;
    }

    this->roots = intervals;
    this->rootsNum = rootNum;
}

Interval NonLinear::bisect(int rootPos, double h) {
    Interval curInterval = this->roots.at(rootPos);

    if (func(curInterval.start) == 0) {
        curInterval.end = curInterval.start;
        return curInterval;
    }
    if (func(curInterval.end) == 0) {
        curInterval.start = curInterval.end;
        return curInterval;
    }
    double x0 = curInterval.start;
    double x = curInterval.end;
    double xi = 0;
    int curRoots = 0;

    while (std::abs(x0 - x) > h) {
        xi = x + (x0 - x) / 2;
        if(func(x) < 0 && func(xi) > 0 || func(x) > 0 && func(xi) < 0) {
            x0 = xi;
            curRoots++;
        } else {
            x = xi;
        }
    }

    curInterval.start = x0;
    curInterval.end = x;
    roots.at(rootPos) = curInterval;

    return curInterval;
}

double NonLinear::newton(int rootPos, double h, double maxIteration) {
    Interval curInterval = this->roots.at(rootPos);
    if (func(curInterval.start) * func(curInterval.end) > 0) {
        return 666;
    }

    double x0 = curInterval.start;
    int iter = 0;
    if (func(x0) * derivative(derivative(x0)) > 0) {
        x0 = curInterval.start;
    } else {
        x0 = curInterval.end;
    }
    double x = x0 - func(x0) / derivative(x0);
    while (std::abs(x0 - x) > h) {
        x0 = x;
        x = x - func(x0) / derivative(x0);
        iter++;
        if (iter > maxIteration) {
            break;
        }
    }

    return x;
}

std::vector<double> NonLinear::solve(double h, double maxIteration) {
    std::vector<double> solutions;

    for (int i = 0; i < roots.size(); i++) {
        bisect(i, 10e-4);
        solutions.emplace_back(newton(i, h, maxIteration));
    }

    return solutions;
}

double NonLinear::derivative(double x) {
    double h = 1e-7;
    return (func(x + h) - func(x - h)) / 2 / h;
}

const Interval &NonLinear::getInterval() const {
    return interval;
}

int NonLinear::getRootsNum() const {
    return rootsNum;
}

const std::vector<Interval> &NonLinear::getRoots() const {
    return roots;
}