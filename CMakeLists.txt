cmake_minimum_required(VERSION 3.27)
project(MagnettoPrjct LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory( vtu11 )

add_executable(MagnettoPrjct
        main.cpp
        LinOp.h
        FileIO.h
        MHDProblem1D.h
        MHDProblem1D.cpp
        MHDSolver1D.cpp
        MHDSolver1D.h
        ShockTube1D.cpp
        ShockTube1D.h
        MHDSolver2D.cpp
        MHDSolver2D.h
        MHDSolver2D.h
        NetGeometry.cpp
        NetGeometry.h
        NamesNConstants.h
        NamesNConstants.h)
target_link_libraries( MagnettoPrjct PRIVATE vtu11::vtu11 )