#include "MHDSolver1D.h"
#include "ShockTube1D.h"
#include "NetGeometry.h"
#include "MHDSolver2D.h"

#include <iostream>
#include <format>
#include <string>
#include <string_view>
#include<random>
#include <cassert>

void tests1D(){
    std::vector<double> BrioWu_L1{1., 0., 0., 0., 1., 0.75, 1., 0.};
    std::vector<double> BrioWu_R1{0.125, 0., 0., 0., 0.1, 0.75, -1., 0.};
    /*rho  u     v    w     p    Bx                             By                        Bz*/
    std::vector<double> BrioWu_L2{1.08, 1.2, 0.01, 0.5, 0.95, 4./(std::sqrt(4.*M_PI)), 3.6/(std::sqrt(4.*M_PI)), 2./(std::sqrt(4.*M_PI))};
    std::vector<double> BrioWu_R2{1., 0., 0., 0., 1., 4./(std::sqrt(4.*M_PI)), 4./(std::sqrt(4.*M_PI)), 2./(std::sqrt(4.*M_PI))};

    //double gam_hcr_1 = 5./3.;
    double gam_hcr_1 = 2.0;
    double x0_1 = -0.5;
    double X_1 = 0.5;
    bool what_is_L_1 = false;
    double t0_1 = 0.;
    //double T_1 = 0.2;
    double T_1 = 0.1;
    double h_1 = 0.0025;
    double tau_1 = 0.00002;
    double gam_courant = 0.1;
    //double gam_courant = 0.2;
    MHDProblem1D problem1(gam_hcr_1, x0_1, X_1, t0_1, T_1, h_1, tau_1, gam_courant, what_is_L_1);

    ShockTube1D BrioWuTest = ShockTube1D(gam_hcr_1, BrioWu_L2, BrioWu_R2, (problem1.X+problem1.x0)/2);

    problem1.initStateFunc = BrioWuTest.initDistrib;

    problem1.leftBoundaryFunction = BrioWuTest.leftBound;

    problem1.rightBoundaryFunction = BrioWuTest.rightBound;

    //HLLScheme(problem1);

    double h2 = 0.0025;
    double tau2 = 0.00002;
    MHDProblem1D problem2(gam_hcr_1, x0_1, X_1, t0_1, T_1, h2, tau2, gam_courant, what_is_L_1);
    problem2.initStateFunc = BrioWuTest.initDistrib;
    problem2.leftBoundaryFunction = BrioWuTest.leftBound;
    problem2.rightBoundaryFunction = BrioWuTest.rightBound;
    //HLLCScheme(problem2);

    double h3 = 0.0025;
    double tau3 = 0.00002;
    MHDProblem1D problem3(gam_hcr_1, x0_1, X_1, t0_1, T_1, h3, tau3, gam_courant, what_is_L_1);
    problem3.initStateFunc = BrioWuTest.initDistrib;
    problem3.leftBoundaryFunction = BrioWuTest.leftBound;
    problem3.rightBoundaryFunction = BrioWuTest.rightBound;
    //HLLDScheme(problem3);
//!!!!!!!
    std::vector<int> space_numbers{64,128,256,512,1024};
    for(int k = 0; k < 5; ++k){
        double hh = (X_1 - x0_1)/space_numbers[k];
        MHDProblem1D problem(gam_hcr_1, x0_1, X_1, t0_1, T_1, hh, tau_1, gam_courant, what_is_L_1);

        ShockTube1D BrioWu = ShockTube1D(gam_hcr_1, BrioWu_L1, BrioWu_R1, (problem.X+problem.x0)/2);

        problem.initStateFunc = BrioWu.initDistrib;

        problem.leftBoundaryFunction = BrioWu.leftBound;

        problem.rightBoundaryFunction = BrioWuTest.rightBound;
        HLLScheme(problem, std::format("BrioNetTest_HLL_{}", k));
        HLLCScheme(problem, std::format("BrioNetTest_HLLC_{}", k));
        HLLDScheme(problem, std::format("BrioNetTest_HLLD_{}", k));
    }

//!!!!!!!!!!!!!!!
//    // Распространение циркулярно-поляризованной альфвеновской волны
//    //double alpha_angle = M_PI/6;
//    double gam_hcr_2 = 5./3.;
//    double x0_2 = 0.;
//   // double X_2 = 1/std::cos(alpha_angle); //~1.1547
//    double X_2 = 1.;
//    bool what_is_L_2 = false;
//    double t0_2 = 0.;
//    double T_2 = 0.2;
//    //double h_2 = 0.0180422;
//    double h_2 = 0.0090211;
//    double tau_2 = 0.002;
//    double gam_courant_2 = 0.85;
//
//    std::vector<int> space_numbers{8,16,32,64,128};
//    for(int k = 0; k < 5; ++k){
//        double hh = (X_2 - x0_2)/space_numbers[k];
//        MHDProblem1D problemAlfven(gam_hcr_2,x0_2, X_2, t0_2, T_2, hh, tau_2, gam_courant_2, what_is_L_2);
//        problemAlfven.initStateFunc =  ([&](double x) {
//            double rho_2 = 1.0;
//            double u_2 = 0.;
//            std::function<double(double)> v_2 = ([&](double x) {return 0.1 * std::sin(2*M_PI*x);});
//            std::function<double(double)> w_2 = ([&](double x) {return 0.1 * std::cos(2*M_PI*x);});
//            double p_2 = 0.1;
//            double Bx_2 = std::sqrt(4*M_PI);
//            std::function<double(double)> By_2 = ([&](double x) {return 0.1 *std::sqrt(4*M_PI)*std::sin(2*M_PI*x);});
//            std::function<double(double)> Bz_2 = ([&](double x) {return 0.1 * std::sqrt(4*M_PI)*std::cos(2*M_PI*x);});
//            return state_from_primitive_vars(rho_2, u_2, v_2(x), w_2(x), p_2, Bx_2, By_2(x), Bz_2(x), gam_hcr_2);
//        });
//        problemAlfven.periodic_boundaries = true;
//        //problemAlfven.leftBoundaryFunction = ([&](double t){return problemAlfven.initStateFunc(problemAlfven.x0);});
//        //problemAlfven.rightBoundaryFunction = ([&](double t){return problemAlfven.initStateFunc(problemAlfven.X);});
//        HLLScheme(problemAlfven, std::format("AlfvenWaveTest_HLL_{}", k));
//        HLLCScheme(problemAlfven, std::format("AlfvenWaveTest_HLLC_{}", k));
//        HLLDScheme(problemAlfven, std::format("AlfvenWaveTest_HLLD_{}", k));
//    }
}


/* ТЕСТЫ ГЕОМЕТРИИ*/
void testNode() {
    Node n(0, 1.0, 2.0, 3.0);
    assert(n.ind == 0);
    assert(n.x == 1.0);
    assert(n.y == 2.0);
    assert(n.z == 3.0);
    std::cout << "Node test passed!" << std::endl;
}
void testElement() {
    std::vector<int> nodeIndexes = {0, 1, 2};
    Element e(0, nodeIndexes, 3);
    assert(e.ind == 0);
    assert(e.nodeIndexes == nodeIndexes);
    assert(e.dim == 3);
    std::cout << "Element test passed!" << std::endl;
}
void testEdge() {
    std::vector<double> normalVec = {0.0, 1.0};
    std::vector<double> midP = {0.5, 1.5};
    Edge e(0, 0, 1, -1, -1, 1, 1.0, normalVec, midP);
    assert(e.ind == 0);
    assert(e.nodeInd1 == 0);
    assert(e.nodeInd2 == 1);
    assert(e.neighbourInd1 == -1);
    assert(e.neighbourInd2 == -1);
    assert(e.orientation == 1);
    assert(e.length == 1.0);
    assert(e.normalVector == normalVec);
    assert(e.midPoint == midP);
    std::cout << "Edge test passed!" << std::endl;
}
void testAreaCalc() {
    std::vector<Node> nodes = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 0.0, 0.0),
            Node(2, 0.0, 1.0, 0.0)
    };
    std::vector<int> nodeIndexes = {0, 1, 2};
    Element e(0, nodeIndexes, 3);
    NodePool np(3, nodes);
    double area = areaCalc(e, np);
    assert(area == 0.5); // Area of a triangle with base and height of 1
    std::cout << "Area calculation test passed!" << std::endl;
}
void testCentroidCalc() {
    std::vector<Node> nodes = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 0.0, 0.0),
            Node(2, 0.0, 1.0, 0.0)
    };
    std::vector<int> nodeIndexes = {0, 1, 2};
    Element e(0, nodeIndexes, 3);
    NodePool np(3, nodes);
    std::vector<double> centroid = getElementCentroid2D(e, np);
    assert(centroid[0] == 1.0 / 3.0);
    assert(centroid[1] == 1.0 / 3.0);
    std::cout << "Centroid calculation test passed!" << std::endl;
}
void testEdgePool() {
    std::vector<Node> nodes = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 0.0, 0.0),
            Node(2, 0.0, 1.0, 0.0),
            Node(3, 1.0, 1.0, 0.0)
    };
    std::vector<Element> elements = {
            Element(0, {0, 1, 2}, 3),
            Element(1, {1, 2, 3}, 3)
    };
    NodePool np(4, nodes);
    ElementPool ep(3, 2, elements);
    EdgePool edgePool(np, ep);
    assert(edgePool.edges.size() == 4);  // Check if the expected number of edges is created
    std::cout << "EdgePool test passed!" << std::endl;
}
void testNeighbourService() {
    std::vector<Node> nodes = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 0.0, 0.0),
            Node(2, 0.0, 1.0, 0.0),
            Node(3, 1.0, 1.0, 0.0)
    };
    std::vector<Element> elements = {
            Element(0, {0, 1, 2}, 3),
            Element(1, {1, 2, 3}, 3)
    };
    NodePool np(4, nodes);
    ElementPool ep(3, 2, elements);
    EdgePool edgePool(np, ep);
    NeighbourService ns(np, ep, edgePool);

    std::unordered_set<int> nodeNeighbours = ns.getNodeNeighbours(1);
    assert(nodeNeighbours.size() == 2);  // Node 1 should be shared by two elements
    std::cout << "NeighbourService test passed!" << std::endl;
}

int main() {
    testNode();
    testElement();
    testEdge();
    testAreaCalc();
    testCentroidCalc();
    testEdgePool();
    testNeighbourService();

    World world("InputData/mesh005.txt");
    //world.display();

    MHDSolver2D solver(world);
    solver.runSolver();
    writeVTU("OutputData/2D/output005.vtu", world, solver.elemUs);
    std::cout << "solver complete" << std::endl;

    return 0;
}
