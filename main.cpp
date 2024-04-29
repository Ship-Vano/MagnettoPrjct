#include <iostream>
#include "MHDSolver.h"

int main() {

    double gam_hcr = 5./3.;
    double x0 = -0.5;
    double X = 0.5;
    bool what_is_L = false;
    double t0 = 0.;
    double T = 0.1;
    double h = 0.0025;
    double tau = 0.00002;

    std::function<std::vector<double>(double)> BrioWu1 = ([&](double x) {
        if(x < 0.){
            const double rhoL = 1.;
            const double uL = 0.;
            const double vL = 0.;
            const double wL = 0.;
            const double pL = 1.;
            const double BxL = 0.75;
            const double ByL = 1.;
            const double BzL = 0.;
            return state_from_primitive_vars(rhoL, uL, vL, wL, pL, BxL, ByL, BzL, gam_hcr);
        }
        else{
            const double rhoR = 0.125;
            const double uR = 0.;
            const double vR = 0.;
            const double wR = 0.;
            const double pR = 0.1;
            const double BxR = 0.75;
            const double ByR = -1.;
            const double BzR = 0.;
            return state_from_primitive_vars(rhoR, uR, vR, wR, pR, BxR, ByR, BzR, gam_hcr);
        }
    });

    std::function<std::vector<double>(double)> leftBound1 =([&](double t){
        const double rhoL = 1.;
        const double uL = 0.;
        const double vL = 0.;
        const double wL = 0.;
        const double pL = 1.;
        const double BxL = 0.75;
        const double ByL = 1.;
        const double BzL = 0.;
        return state_from_primitive_vars(rhoL, uL, vL, wL, pL, BxL, ByL, BzL, gam_hcr);
    });

    std::function<std::vector<double>(double)> rightBound1 = ([&](double t){
        const double rhoR = 0.125;
        const double uR = 0.;
        const double vR = 0.;
        const double wR = 0.;
        const double pR = 0.1;
        const double BxR = 0.75;
        const double ByR = -1.;
        const double BzR = 0.;
        return state_from_primitive_vars(rhoR, uR, vR, wR, pR, BxR, ByR, BzR, gam_hcr);
    });

    MHDProblem problem1(gam_hcr, x0, X, t0, T, h, tau, what_is_L);

    problem1.initStateFunc = BrioWu1;

    problem1.leftBoundaryFunction = leftBound1;

    problem1.rightBoundaryFunction = rightBound1;

    //HLLScheme(problem1);
    double h2 = 0.005;
    double tau2 = 0.0002;
    MHDProblem problem2(gam_hcr, x0, X, t0, T, h2, tau2, what_is_L);
    problem2.initStateFunc = BrioWu1;
    problem2.initStateFunc = BrioWu1;
    problem2.leftBoundaryFunction = leftBound1;
    problem2.rightBoundaryFunction = rightBound1;
    HLLCScheme(problem2);
    return 0;
}
