cmake_minimum_required(VERSION 3.27)
project(MagnettoPrjct)

set(CMAKE_CXX_STANDARD 20)

add_executable(MagnettoPrjct
        main.cpp
        LinOp.h
        FileIO.h
        MHDProblem.h
        MHDProblem.cpp
        MHDSolver.cpp
        MHDSolver.h
        ShockTube1D.cpp
        ShockTube1D.h)
