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

    MHDProblem problem1(gam_hcr, x0, X, t0, T, h, tau, what_is_L);

    problem1.initStateFunc = ([&](double x) {
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

    problem1.leftBoundaryFunction = ([&](double t){
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

    problem1.rightBoundaryFunction = ([&](double t){
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

    HLLScheme(problem1);
    return 0;
}
