cmake_minimum_required(VERSION 3.16)
project(mv_02)

set(CMAKE_CXX_STANDARD 20)

add_executable(mv_02 main.cpp matrix.cpp matrix.h vector.cpp vector.h nonLinear.cpp nonLinear.h eigen.cpp eigen.h)
