//
// Created by Иван on 10/21/2024.
//

#ifndef MAGNETTOPRJCT_MHDSOLVER2D_H
#define MAGNETTOPRJCT_MHDSOLVER2D_H

#include "NetGeometry.h"
#include "MHDSolver1D.h"

class MHDSolver2D {
public:
    // mesh
    World geometryWorld;

    // physical quantities
    double gam_hcr = 2.0;
    double startTime = 0.0;  // время отсчёта
    double finalTime = 0.1;   // время окончания
    double tau = 0.0001; // шаг по времени
    double cflNum = 0.4; // число Куранта
    int MAX_ITERATIONS = 10000;

    //states
    std::vector<std::vector<double>> nodeUs; // state U in nodes
    std::vector<std::vector<double>> elemUs; // state U in elements
    std::vector<std::vector<double>> edgeUs; // state U in edges
    std::vector<std::vector<double>> initElemUs;
    std::vector<double> initBns;
    std::vector<double> bNs; //Bns at edges

    static std::vector<double> rotateStateFromNormalToAxisX(std::vector<double>& U, const std::vector<double>& n);
    static std::vector<double> rotateStateFromAxisToNormal(std::vector<double>& U, const std::vector<double>& n);
    void setInitElemUs();
    void runSolver();

    // Начальное состояние системы
    std::function<std::vector<double>(double)> initStateFunc;
    bool initStateFunc_is_set = false;

    // Левые граничные условия
    std::function<std::vector<double>(double)> leftBoundaryFunction;
    bool leftBoundaryFunction_is_set = false;

    // Правые граничные условия
    std::function<std::vector<double>(double)> rightBoundaryFunction;
    bool rightBoundaryFunction_is_set = false;

    MHDSolver2D(const World& world);
};

void writeVTU(const std::string& filename, const World& geometryWorld, const std::vector<std::vector<double>>& elemUs);


/*                       0      1      2      3    4   5   6   7
 * U (general state):  rho,  rho*u, rho*v, rho*w,  e,  Bx, Bz, By
 * gasU                rho,  rho*u, rho*v, rho*w,  e
 * magU                 Bx,    Bz,    By
 * */

//void solverHLL2D(const World& world);

#endif //MAGNETTOPRJCT_MHDSOLVER2D_H
