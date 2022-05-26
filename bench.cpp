#include <iostream>
#include <cmath>
#include <random>
#include <complex>
#include <chrono>

#include "matrix.h"
#include "eigen.h"

#include "Eigen/Dense"

using Eigen::MatrixXd;

void benchTime(int minSize, int maxSize, int sizeStep, double randMin, double randMax, int repeatNum) {

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int size = minSize; size <= maxSize; size += sizeStep) {
        long long timer = 0;
        int invalid = 0;
        int valid = 0;
        int iterations = 0;

        for (int repeat = 0; repeat < repeatNum; repeat++) {
            std::uniform_real_distribution<double> distr(randMin, randMax);

            Matrix matrix = Matrix(70);
            MatrixXd control(70, 70);
            for (int i = 0; i < 70; i++) {
                for (int j = 0; j < 70; j++) {
                    double randNum = distr(gen);
                    matrix(i, j) = randNum;
                    control(i, j) = randNum;
                }
            }

            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            eigen eigenvv = eigen::powerIteration(matrix, Vector(size, 1), 1e-8, size);
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            timer += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            Eigen::ComplexEigenSolver<Eigen::MatrixXd> eigenSolver(control);
            if (eigenSolver.info() != Eigen::Success) {
                std::cout << "EIGENSOLVER DID NOT FOUND VECTORS" << std::endl;
            }
            else if (eigenvv.getEigenKind() == EigenKind::NoValues) {
                invalid++;
                //timer -= std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            }
            else {

                double diff = std::abs(eigenSolver.eigenvalues()[eigenSolver.eigenvalues().size() - 1])
                              - std::abs(eigenvv.getEigenValues()[eigenvv.getNumberOfValues() - 1]);
                if (std::abs(diff) > 1e-6) {
                    std::cout << "EIGEN: " << eigenSolver.eigenvalues()[eigenSolver.eigenvalues().size() - 1];
                    std::cout << " OUR: " << eigenvv.getEigenValues()[eigenvv.getNumberOfValues() - 1];
                    std::cout << " DIFF: " << diff << std::endl;
                    invalid++;
                } else {
                    valid++;
                }
            }

            iterations += eigenvv.getIterations();
        }

        std::cout << size << " " << timer / repeatNum << " " << (double)invalid / repeatNum * 100 << " " << iterations / repeatNum << std::endl;
    }
}

void benchQR(int minSize, int maxSize, int sizeStep, double randMin, double randMax, int repeatNum) {

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int size = minSize; size <= maxSize; size += sizeStep) {
        long long timer = 0;
        int invalid = 0;
        int valid = 0;
        int iterations = 0;

        for (int repeat = 0; repeat < repeatNum; repeat++) {
            std::uniform_real_distribution<double> distr(randMin, randMax);

            Matrix matrix = Matrix(size);
            MatrixXd control(size, size);
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    double randNum = distr(gen);
                    matrix(i, j) = randNum;
                    control(i, j) = randNum;
                }
            }

            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            eigen eigenvv = eigen::QR(matrix);
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            timer += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            if (eigenvv.getEigenKind() == EigenKind::NoValues || eigenvv.getIterations() == 20000) {
                invalid++;
                timer -= std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            }

            iterations += eigenvv.getIterations();
        }

        std::cout << size << " " << timer / repeatNum << " " << (double)invalid / repeatNum * 100 << " " << iterations / repeatNum << std::endl;

    }
}

void benchStability(int size, double randMin, double randMax, double randStep, int repeatNum) {

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int r = randMin; r <= randMax; r += randStep) {
        long long timer = 0;
        int invalid = 0;
        int valid = 0;

        for (int repeat = 0; repeat < repeatNum; repeat++) {
            std::uniform_real_distribution<double> distr(r - randStep, r);

            Matrix matrix = Matrix(size);
            MatrixXd control(size, size);
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    double randNum = distr(gen);
                    matrix(i, j) = randNum;
                    control(i, j) = randNum;
                }
            }

            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            eigen eigenvv = eigen::powerIteration(matrix, Vector(size, 1), 1e-8, 100000);
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            timer += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            Eigen::ComplexEigenSolver<Eigen::MatrixXd> eigenSolver(control);
            if (eigenSolver.info() != Eigen::Success) {
                std::cout << "EIGENSOLVER DID NOT FOUND VECTORS" << std::endl;
            }
            else if (eigenvv.getEigenKind() == EigenKind::NoValues) {
                invalid++;
                //timer -= std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            }
            else {

                double diff = std::abs(eigenSolver.eigenvalues()[eigenSolver.eigenvalues().size() - 1])
                              - std::abs(eigenvv.getEigenValues()[eigenvv.getNumberOfValues() - 1]);
                if (std::abs(diff) > 1e-6) {
                    std::cout << "EIGEN: " << eigenSolver.eigenvalues()[eigenSolver.eigenvalues().size() - 1];
                    std::cout << " OUR: " << eigenvv.getEigenValues()[eigenvv.getNumberOfValues() - 1];
                    std::cout << " DIFF: " << diff << std::endl;
                    invalid++;
                } else {
                    valid++;
                }
            }
        }

        std::cout << r << " " << timer / repeatNum << " " << (double)invalid / repeatNum * 100 << std::endl;

    }
}