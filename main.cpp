#include <iostream>
#include "MHDSolver.h"
#include "ShockTube1D.h"

int main() {

    std::vector<double> BrioWu_L1{1., 0., 0., 0., 1., 0.75, 1., 0.};
    std::vector<double> BrioWu_R1{0.125, 0., 0., 0., 0.1, 0.75, -1., 0.};
                                  /*rho  u     v    w     p    Bx                             By                        Bz*/
    std::vector<double> BrioWu_L2{1.08, 1.2, 0.01, 0.5, 0.95, 4./(std::sqrt(4.*M_PI)), 3.6/(std::sqrt(4.*M_PI)), 2./(std::sqrt(4.*M_PI))};
    std::vector<double> BrioWu_R2{1., 0., 0., 0., 1., 4./(std::sqrt(4.*M_PI)), 4./(std::sqrt(4.*M_PI)), 2./(std::sqrt(4.*M_PI))};

    double gam_hcr = 5./3.;
    double x0 = -0.5;
    double X = 0.5;
    bool what_is_L = false;
    double t0 = 0.;
    double T = 0.2;
    double h = 0.0025;
    double tau = 0.00002;
    MHDProblem problem1(gam_hcr, x0, X, t0, T, h, tau, what_is_L);


    ShockTube1D BrioWuTest = ShockTube1D(gam_hcr, BrioWu_L2, BrioWu_R2, (problem1.X+problem1.x0)/2);


    problem1.initStateFunc = BrioWuTest.initDistrib;

    problem1.leftBoundaryFunction = BrioWuTest.leftBound;

    problem1.rightBoundaryFunction = BrioWuTest.rightBound;

    HLLScheme(problem1);

    double h2 = 0.0025;
    double tau2 = 0.00002;
    MHDProblem problem2(gam_hcr, x0, X, t0, T, h2, tau2, what_is_L);
    problem2.initStateFunc = BrioWuTest.initDistrib;
    problem2.leftBoundaryFunction = BrioWuTest.leftBound;
    problem2.rightBoundaryFunction = BrioWuTest.rightBound;
    HLLCScheme(problem2);

    double h3 = 0.0025;
    double tau3 = 0.00002;
    MHDProblem problem3(gam_hcr, x0, X, t0, T, h3, tau3, what_is_L);
    problem3.initStateFunc = BrioWuTest.initDistrib;
    problem3.leftBoundaryFunction = BrioWuTest.leftBound;
    problem3.rightBoundaryFunction = BrioWuTest.rightBound;
    HLLDScheme(problem3);
    return 0;
}
