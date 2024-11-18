//
// Created by Иван on 10/21/2024.
//

#ifndef MAGNETTOPRJCT_MHDSOLVER2D_H
#define MAGNETTOPRJCT_MHDSOLVER2D_H

#include "FileIO.h"
#include "NetGeometry.h"

class MHDSolver2D {

    // mesh
    World geometryWorld;

    // physical quantities
    std::vector<std::vector<double>> nodeGasUs; // gas U in nodes
    std::vector<std::vector<double>> elemGasUs; // gas U in elements
    std::vector<std::vector<double>> edgeGasUs; // gas U in edges
    std::vector<std::vector<double>> nodeMagUs; // magnetic U in nodes
    std::vector<std::vector<double>> elemMagUs; // magnetic U in elements
    std::vector<std::vector<double>> edgeMagUs; // magnetic U in edges

    void runSolver();
};



/*                       0      1      2      3    4   5   6   7
 * U (general state):  rho,  rho*u, rho*v, rho*w,  e,  Bx, Bz, By
 * gasU                rho,  rho*u, rho*v, rho*w,  e
 * magU                 Bx,    Bz,    By
 * */

//void solverHLL2D(const World& world);

#endif //MAGNETTOPRJCT_MHDSOLVER2D_H
