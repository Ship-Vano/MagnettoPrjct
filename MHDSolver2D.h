//
// Created by Иван on 10/21/2024.
//

#ifndef MAGNETTOPRJCT_MHDSOLVER2D_H
#define MAGNETTOPRJCT_MHDSOLVER2D_H

#include "NetGeometry.h"
#include "MHDSolver1D.h"

class MHDSolver2D {

    // mesh
    World geometryWorld;

    // physical quantities
    double gam_hcr = 3.0/2.0;
    //states
    std::vector<std::vector<double>> nodeUs; // state U in nodes
    std::vector<std::vector<double>> elemUs; // state U in elements
    std::vector<std::vector<double>> edgeUs; // state U in edges

    //fluxes
    std::vector<std::vector<double>> fluxes;    // MHD (HLLD) fluxes (from one element to another "<| -> |>")

    std::vector<double> rotateStateFromNormalToAxisX(std::vector<double>& U, const std::vector<double>& n);
    std::vector<double> rotateStateFromAxisToNormal(std::vector<double>& U, const std::vector<double>& n);
    void runSolver();
};



/*                       0      1      2      3    4   5   6   7
 * U (general state):  rho,  rho*u, rho*v, rho*w,  e,  Bx, Bz, By
 * gasU                rho,  rho*u, rho*v, rho*w,  e
 * magU                 Bx,    Bz,    By
 * */

//void solverHLL2D(const World& world);

#endif //MAGNETTOPRJCT_MHDSOLVER2D_H
