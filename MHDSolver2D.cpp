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


/*вращения*/
// n = {cos(j), sin(j)} = {n.x, n.y}
// ({cos(j), -sin(j), 0}, {sin(j), cos(j), 0}, {0,0,1}) --- around OZ counterclockwise
// rotate 1: from normal to OX:      {v.x * n.x + vy * n.y, - v.x * n.y + v.y * n.x, v.z}
// вращаем {u,v,w} и {Bx, By, Bz}, остальные остаются на месте
std::vector<double> MHDSolver2D::rotateStateFromAxisToNormal(vector<double> &U, const vector<double>& n) {
    std::vector<double> res(U);
    res[1] =  U[1]*n[0] + U[2]*n[1];
    res[2] = -U[1]*n[1] + U[2]*n[0];
    res[5] =  U[5]*n[0] + U[6]*n[1];
    res[6] = -U[5]*n[1] + U[6]*n[0];
    return res;
}

// ({cos(j), sin(j), 0}, {-sin(j), cos(j), 0}, {0,0,1}) --- around OZ clockwise
// rotate 2: from OX to normal:      {v.x * n.x - vy * n.y,  v.x * n.y + v.y * n.x, v.z}
std::vector<double> MHDSolver2D::rotateStateFromNormalToAxisX(vector<double> &U, const vector<double>& n) {
    std::vector<double> res(U);
    res[1] =  U[1]*n[0] - U[2]*n[1];
    res[2] =  U[1]*n[1] + U[2]*n[0];
    res[5] =  U[5]*n[0] - U[6]*n[1];
    res[6] =  U[5]*n[1] + U[6]*n[0];
    return res;
}

void MHDSolver2D::runSolver() {
    EdgePool edgePool = geometryWorld.getEdgePool();
    ElementPool elPool = geometryWorld.getElementPool();
    NodePool nodePool = geometryWorld.getNodePool();
    NeighbourService ns = geometryWorld.getNeighbourService();

    // сделать старые дубликаты состояний (предыдущие состояния) чтобы в новые записывать расчёты
    std::vector<std::vector<double>> elemUs_prev(elemUs);
    std::vector<std::vector<double>> nodeUs_prev(nodeUs);
    std::vector<std::vector<double>> edgeUs_prev(edgeUs);

    std::vector<std::vector<double>> fluxes_temp(edgePool.edgeCount, std::vector<double>(8,0.0));


    double tau = 0.0001; //TODO: подбирать tau из условия Куранта

    //(2) вычисляем потоки, проходящие через каждое ребро
    int count = 0; // TODO: для использования параллелек нужно в конструкторе проинициализировать флаксы сразу и разобраться, как брать номер итератора. мб вообще убрать итератор и стандартно по i делать цикл по эджам))
    for(const auto& edge : edgePool.edges){
        Node node1 = nodePool.getNode(edge.nodeInd1);
        Node node2 = nodePool.getNode(edge.nodeInd2);
        int neighbour1 = edge.neighbourInd1;
        int neighbour2 = edge.neighbourInd2;
        std::vector<double> U1 = rotateStateFromNormalToAxisX(elemUs[neighbour1], edge.normalVector);
        if(neighbour2 > 0) {
            std::vector<double> U2 = rotateStateFromNormalToAxisX(elemUs[neighbour2], edge.normalVector);
            //fluxes.emplace_back(HLLD_flux(U1, U2, gam_hcr));
            fluxes_temp[count] = HLLD_flux(U1, U2, gam_hcr);
            rotateStateFromAxisToNormal(fluxes_temp[count] , edge.normalVector);
        }
        else{ //здесь задавать гран условие ещё мб на поток
            fluxes_temp[count] = HLLD_flux(U1, U1, gam_hcr);
            rotateStateFromAxisToNormal(fluxes_temp[count] , edge.normalVector);
        }
        ++count;
    }

    //по явной схеме обновляем газовые величины
    for(const auto& elem: elPool.elements){
        int i = elem.ind;
        std::vector<double> fluxSum(8, 0.0);
        for(int edgeIndex : elem.edgeIndexes){
            Edge edge_j = edgePool.edges[edgeIndex];
            if(edge_j.neighbourInd1 == i){
                fluxSum = fluxSum + edge_j.length * fluxes_temp[edgeIndex];
            }
            else{
                fluxSum = fluxSum - edge_j.length * fluxes_temp[edgeIndex];
            }
        }
        elemUs[i] = elemUs_prev[i] - tau/elem.area * fluxSum;
    }

    //корректируем магитные величины
    //находим узловые значения нужных магнитных разностей //(v x B)z в узлах
    std::vector<double> nodeMagDiffs(nodePool.nodeCount, 0.0); //(v x B)z в узлах
    for(const auto& node: nodePool.nodes){
        int tmp_count = 0;
        for(const auto& neighbourEdgeInd: ns.getEdgeNeighborsOfNode(node.ind)){
            nodeMagDiffs[node.ind] += fluxes_temp[neighbourEdgeInd][6];
            ++tmp_count;
        }
        if(tmp_count){
            nodeMagDiffs[node.ind] /= tmp_count;
        }
    }

    //находим новое значение Bn в ребре
    std::vector<double> bNs_prev(bNs);
    for(int i = 0; i < edgePool.edgeCount; ++i){
        Edge edge = edgePool.edges[i];
        bNs[i] = bNs_prev[i] - tau / edge.length * (nodeMagDiffs[edge.nodeInd2] - nodeMagDiffs[edge.nodeInd1]); //возможно, здесь знак надо учесть
    }


    //сносим Bn в центр элемента
    for(const auto& elem: elPool.elements){
        elemUs[elem.ind][5] = 0.0;
        elemUs[elem.ind][6] = 0.0;
        std::vector<double> centroid = getElementCentroid2D(elem, nodePool);
        for(const auto& edgeInd: elem.edgeIndexes){
            Edge edge = edgePool.edges[edgeInd];
            if(edge.neighbourInd1 == elem.ind) {
                // у первого соседа в эдже заданы ноды в порядке полодительного обхода и нормаль тоже
                const auto nodeInElemeInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(), edge.nodeInd1);
                int node_before_ind =  nodeInElemeInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim-1] :  *(nodeInElemeInd-1);
                //        auto itNode1 = std::find(elementNodes.begin(), elementNodes.end(), node1);
                Node node_before = nodePool.getNode(node_before_ind);
                elemUs[elem.ind][5] += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                elemUs[elem.ind][6] += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
            }
            else{
                // а вот для второго нужно умножать на -1 и в обратном порядке
                const auto nodeInElemeInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(), edge.nodeInd2);
                int node_before_ind =  nodeInElemeInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim-1] :  *(nodeInElemeInd-1);
                //        auto itNode1 = std::find(elementNodes.begin(), elementNodes.end(), node2);
                Node node_before = nodePool.getNode(node_before_ind);
                elemUs[elem.ind][5] -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                elemUs[elem.ind][6] -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
            }
        }
    }
}

// (1)подбираем оптимальное число шага по времени (Курант хард ту килл)
// (2)вычисляем потоки все на рёбрах
// (3)далее идём по всем элементам
// (3.1)вращаем газ (онли скорости), вращаем магнит
// (3.2)собираем два состояния и кидаем в HLLD
// (3.3)полученный поток вращаем обратно, записываем в потоки на рёбрах (проекция потока на элементе на нормаль одного и ребра это поток на ребре, умноженный на ориентацию)
// (3.4)дальше явная схема: д/дт (u_i) * s_i + Sum(j = 1, 3) l_j vec{n}_j F_i = 0
// (3.5)вычислили новое состояние, далее газ не трогаем
// (3.6)теперь апдейтим магнит (чтобы чистенкая дивергенция была):
// (3.6.1)д/дт (B_n) = (v_x*B_y - v_y*B_x)_a - (v_x*B_y - v_y*B_x)_b
// здесь компоненты магнита берём из поток на ребрах, которые прилегают к данной ноде
// b_n = b_x * n.x + b_y * n.y
// a - левый узел, b - правый узел ребра элемента
// обновлённый магнит записываем в массивы переменных
// идём дальше по времени (повтор пред шагов)

