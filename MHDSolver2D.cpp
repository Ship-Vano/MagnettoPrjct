//
// Created by Иван on 10/21/2024.
//

#include "MHDSolver2D.h"






//void solverHLL2D(const World& world){
//    NodePool np = world.getNodePool();
//    ElementPool ep = world.getElementPool();
//    std::vector<std::vector<double>> states(np.nodeCount, std::vector<double>(5, 0.0));
//    /*
//     *     states[i] <---> np.nodes[i]
//     *
//     * */
//
//    //test1:   [0, 1] x [0, 1]
//    /*   quadrant            rho        u         v         p
//     *   x>0.5, y>0.5:          0.5313     0.0       0.0       0.4
//     *   x<0.5, y>0.5:          1.0        0.7276    0.0       1.0
//     *   x<0.5, y<0.5:          0.8        0.0       0.0       1.0
//     *   x>0.5, y<0.5:          1.0        0.0       0.7276    1.0
//     *   FREE-FLOW boundary conditions
//     * */
//    //printf("elCount = %d", ep.elCount);
//
//    //initializer
//    double gam_hcr = 3.0/2.0;
//    double w = 0.0;
//    double rho1 = 0.5313; double rho2 = 1.0; double rho3 = 0.8; double rho4 = 1.0;
//    double u1 = 0.0; double u2 = 0.7276; double u3 = 0.0; double u4 = 0.0;
//    double v1 = 0.0; double v2 = 0.0; double v3 = 0.0; double v4 = 0.7276;
//    double p1 = 0.4; double p2 = 1.0; double p3 = 1.0; double p4 = 1.0;
//    for(int i = 0; i < np.nodeCount; ++i){
//        Node tmp_node = np.nodes[i];
//        if(tmp_node.x > 0.5 && tmp_node.y > 0.5){
//           // states[i] = state_from_primitive_vars(rho1, u1, v1, w, p1, gam_hcr);
//        }
//        else if(tmp_node.x < 0.5 && tmp_node.y > 0.5){
//           // states[i] = state_from_primitive_vars(rho2, u2, v2, w, p2, gam_hcr);
//        }
//        else if(tmp_node.x < 0.5 && tmp_node.y < 0.5){
//           // states[i] = state_from_primitive_vars(rho3, u3, v3, w, p3, gam_hcr);
//        }
//        else{
//           // states[i] = state_from_primitive_vars(rho4, u4, v4, w, p4, gam_hcr);
//        }
//    }
//
//
//    //solver
//    for(int i = 0; i < ep.elCount; ++i){
//        //пока делаем для n=4, для n=3 попозже...
//        //std::vector<int> nodeInds = ep.elements[i].nodeIndexes;
//
//    }
//
//}
void MHDSolver2D::runSolver() {
    EdgePool edgePool = geometryWorld.getEdgePool();
    ElementPool elPool = geometryWorld.getElementPool();
    NodePool nodePool = geometryWorld.getNodePool();
    NeighbourService ns = geometryWorld.getNeighbourService();

    for(const auto& edge : edgePool.edges){
        Node node1 = nodePool.getNode(edge.nodeInd1);
        Node node2 = nodePool.getNode(edge.nodeInd2);
        int neighbour1 = edge.neighbourInd1;
        int neighbour2 = edge.neighbourInd2;
        std::vector<double> gasU1 = elemGasUs[neighbour1];
        std::vector<double> gasU2 = elemGasUs[neighbour2];
        std::vector<double> magU1 = elemMagUs[neighbour1];
        std::vector<double> magU2 = elemMagUs[neighbour2];

        // n = {cos(j), sin(j)} = {n.x, n.y}
        // ({cos(j), -sin(j), 0}, {sin(j), cos(j), 0}, {0,0,1}) --- around OZ counterclockwise
        // rotate 1: from normal to OX:      {v.x * n.x + vy * n.y, - v.x * n.y + v.y * n.x, v.z}
        // ({cos(j), sin(j), 0}, {-sin(j), cos(j), 0}, {0,0,1}) --- around OZ clockwise
        // rotate 2: from OX to normal:      {v.x * n.x - vy * n.y,  v.x * n.y + v.y * n.x, v.z}

        // вращаем газ (онли скорости), вращаем магнит
        // собираем два состояния и кидаем в HLLD
        // полученный поток вращаем обратно, записываем в потоки на рёбрах (поток на элементе это поток на ребре, умноженный на ориентацию)
        // дальше явная схема: д/дт (u_i) * s_i + Sum(j = 1, 3) j_j vec{n}_j F_i = 0
        // вычислили новое состояние, далее газ не трогаем
        // теперь апдейтим магнит (чтобы чистенкая дивергенция была):
            // д/дт (B_n) = (v_x*B_y - v_y*B_x)_a - (v_x*B_y - v_y*B_x)_b
            // здесь компоненты магнита берём из узлов (среднее из соседей)
            // b_n = b_x * n.x + b_y * n.y
            // a - левый узел, b - правый узел ребра элемента
    }

}
