#include <iostream>
#include <cmath>
#include <numbers>

double f(double x) {
    return (std::pow(x, 9) + std::numbers::pi)
        * cos(log(std::pow(x, 2) + 1)) /
        std::pow(std::numbers::e, std::pow(x, 2)) - x / 2022;
}

double df(double x) {
    return -2 * std::pow(std::numbers::e, -std::pow(x, 2)) * (std::pow(x, 9) + std::numbers::pi)
           * x * sin(log(std::pow(x, 2) + 1)) / (std::pow(x, 2) + 1)

           -2 * std::pow(std::numbers::e, -std::pow(x, 2)) * (std::pow(x, 9) + std::numbers::pi) * x
            * cos(log(std::pow(x, 2) + 1))

           +9 * std::pow(std::numbers::e, -std::pow(x, 2))
              * std::pow(x, 8) * cos(log(std::pow(x, 2) + 1)) - 1.0/2022;

}

double derivative(double x) {
    double h = 1e-7;
    return (f(x + h) - f(x - h)) / 2 / h /* 1e14*/;
}

struct Interval {
    double start;
    double end;
};

Interval* getNumOfRoots(double a, double b) {
    auto* intervals = new Interval[10]{};
    //double a = -100;
    //double b = 100;
    double h = 0.1;
    int roots = 0;

    double x1 = a;
    double x2 = x1 + h;
    double y1 = f(x1);

    while (x2 < b) {
        double y2 = f(x2);
        if (y1 * y2 < 0) {
            intervals[roots].start = x1;
            intervals[roots].end = x2;
            //std::cout << x1 << " " << x2 << std::endl;
            roots++;
        }

        x1 = x2;
        x2 = x1 + h;
        y1 = y2;
    }

    return intervals;
}

Interval bisect(double a, double b, double h) {
    Interval interval{a, b};
    if (f(a) == 0) {
        interval.end = a;
        return interval;
    }
    if (f(b) == 0) {
        interval.start = b;
        return interval;
    }
    double x0 = a;
    double x = b;
    double xi = 0;
    int roots = 0;

    while (std::abs(x0 - x) > h) {
        xi = x + (x0 - x) / 2;
        if(f(x) < 0 && f(xi) > 0 || f(x) > 0 && f(xi) < 0) {
            x0 = xi;
            roots++;
        } else {
            x = xi;
        }
    }

    interval.start = x0;
    interval.end = x;

    return interval;
}

double newton(double a, double b, double h) {
    if (f(a) * f(b) > 0) {
        return 666;
    }

    double x0 = a;
    int iter = 0;
    if (f(x0) * derivative(derivative(x0)) > 0) {
        x0 = a;
    } else {
        x0 = b;
    }
    double x = x0 - f(x0) / derivative(x0);
    while (std::abs(x0 - x) > h) {
        x0 = x;
        x = x - f(x0) / derivative(x0);
        iter++;
        if (iter > 10000) {
            break;
        }
    }

    /*for(int i = 0; i < 10; i++) {
        x = x - f(x) / derivative(x);
    }*/
    //std::cout << iter << std::endl;

    return x;
}

int main() {
    //std::cout << f(-1.5) << std::endl; // -1.42
    //std::cout << df(-1.5) << std::endl; // 1.85
    //std::cout << derivative(-1.5) << std::endl;

    Interval* intervals = getNumOfRoots(-10, 10);
    for (int i = 0; i < 3; i++) {
        Interval interval = bisect(intervals[i].start, intervals[i].end, 10e-4);
        std::cout << "Using interval [" << interval.start << "; " << interval.end << "]" << std::endl;
        std::cout << newton(interval.start, interval.end, 10e-12) << std::endl;
    }

    return 0;
}
