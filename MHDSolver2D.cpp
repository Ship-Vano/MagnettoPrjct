//
// Created by Иван on 10/21/2024.
//

#include "MHDSolver2D.h"

MHDSolver2D::MHDSolver2D(const World &world): geometryWorld(world), nodeUs(),
    elemUs(), edgeUs(), initElemUs(), initBns(), bNs(){}

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

void MHDSolver2D::setInitElemUs() {
    /*rho  u   v   w   p   Bx   By   Bz gam_hcr*/
    std::vector<double> BrioWu_L1{1.,   0., 0., 0., 1., 0.75, 1., 0., gam_hcr};
    std::vector<double> BrioWu_R1{0.125, 0., 0., 0., 0.1, 0.75, -1., 0., gam_hcr};
    ElementPool ep = geometryWorld.getElementPool();
    initElemUs.resize(ep.elCount, std::vector<double>(8, 0.0));
    for(const auto& elem: ep.elements){
        std::vector<double> centroid = elem.centroid2D;
        int elInd = elem.ind;
        if(centroid[0] < 0.5){
            initElemUs[elInd] = state_from_primitive_vars(BrioWu_L1);
        }
        else {
            initElemUs[elInd] = state_from_primitive_vars(BrioWu_R1);
        }
    }
    EdgePool edgp = geometryWorld.getEdgePool();
    initBns.resize(edgp.edgeCount, 0.0);
    for(const auto& edge: edgp.edges){
        std::vector<double> midP = edge.midPoint;
        int edgeInd = edge.ind;
        std::vector<double> normal = edge.normalVector;
        std::vector<double> state;
        if(midP[0] < 0.5){
            state = state_from_primitive_vars(BrioWu_L1);
        }
        else{
            state = state_from_primitive_vars(BrioWu_R1);
        }
        double Bn = state[5] * normal[0] + state[6] * normal[1];
        initBns[edgeInd] = Bn;
    }
}

void MHDSolver2D::runSolver() {

    setInitElemUs();

    // service
    EdgePool edgePool = geometryWorld.getEdgePool();
    ElementPool elPool = geometryWorld.getElementPool();
    NodePool nodePool = geometryWorld.getNodePool();
    NeighbourService ns = geometryWorld.getNeighbourService();

    // инициализируем состояния
    elemUs = initElemUs;
    bNs = initBns;

    // сделать старые дубликаты состояний (предыдущие состояния) чтобы в новые записывать расчёты
    std::vector<std::vector<double>> elemUs_prev(elemUs);
    std::vector<std::vector<double>> nodeUs_prev(nodeUs);
    std::vector<std::vector<double>> edgeUs_prev(edgeUs);

    // инициализируем вектор потоков через рёбра // MHD (HLLD) fluxes (from one element to another "<| -> |>")
    std::vector<std::vector<double>> fluxes(edgePool.edgeCount, std::vector<double>(8, 0.0));

    double h = edgePool.edges[0].length;

    double currentTime = startTime;

    int iterations = 0;

    while(currentTime < finalTime) {
        elemUs_prev.swap(elemUs);
        nodeUs_prev.swap(nodeUs);
        edgeUs_prev.swap(edgeUs);

        tau = tau_from_cfl(cflNum, h, elemUs, elPool.elCount, gam_hcr);
        currentTime += tau;

        //(2) вычисляем потоки, проходящие через каждое ребро
        int count = 0; // TODO: для использования параллелек нужно в конструкторе проинициализировать флаксы сразу и разобраться, как брать номер итератора. мб вообще убрать итератор и стандартно по i делать цикл по эджам))
        for (const auto &edge: edgePool.edges) {
            Node node1 = nodePool.getNode(edge.nodeInd1);
            Node node2 = nodePool.getNode(edge.nodeInd2);
            int neighbour1 = edge.neighbourInd1;
            int neighbour2 = edge.neighbourInd2;
            std::vector<double> U1 = rotateStateFromNormalToAxisX(elemUs[neighbour1], edge.normalVector);
            if (neighbour2 > 0) {
                std::vector<double> U2 = rotateStateFromNormalToAxisX(elemUs[neighbour2], edge.normalVector);
                //fluxes.emplace_back(HLLD_flux(U1, U2, gam_hcr));
                fluxes[count] = HLLD_flux(U1, U2, gam_hcr);
                rotateStateFromAxisToNormal(fluxes[count], edge.normalVector);
            } else { //здесь задавать гран условие ещё мб на поток
                fluxes[count] = HLLD_flux(U1, U1, gam_hcr);
                rotateStateFromAxisToNormal(fluxes[count], edge.normalVector);
            }
            ++count;
        }

        //по явной схеме обновляем газовые величины
        for (const auto &elem: elPool.elements) {
            int i = elem.ind;
            std::vector<double> fluxSum(8, 0.0);
            for (int edgeIndex: elem.edgeIndexes) {
                Edge edge_j = edgePool.edges[edgeIndex];
                if (edge_j.neighbourInd1 == i) {
                    fluxSum = fluxSum + edge_j.length * fluxes[edgeIndex];
                } else {
                    fluxSum = fluxSum - edge_j.length * fluxes[edgeIndex];
                }
            }
            elemUs[i] = elemUs_prev[i] - tau / elem.area * fluxSum;
        }

        //корректируем магитные величины
        //находим узловые значения нужных магнитных разностей //(v x B)z в узлах
        std::vector<double> nodeMagDiffs(nodePool.nodeCount, 0.0); //(v x B)z в узлах
        for (const auto &node: nodePool.nodes) {
            int tmp_count = 0;
            for (const auto &neighbourEdgeInd: ns.getEdgeNeighborsOfNode(node.ind)) {
                nodeMagDiffs[node.ind] += fluxes[neighbourEdgeInd][6];
                ++tmp_count;
            }
            if (tmp_count) {
                nodeMagDiffs[node.ind] /= tmp_count;
            }
        }

        //находим новое значение Bn в ребре
        std::vector<double> bNs_prev(bNs);
        for (int i = 0; i < edgePool.edgeCount; ++i) {
            Edge edge = edgePool.edges[i];
            bNs[i] = bNs_prev[i] - tau / edge.length * (nodeMagDiffs[edge.nodeInd2] -
                                                        nodeMagDiffs[edge.nodeInd1]); //возможно, здесь знак надо учесть
        }


        //сносим Bn в центр элемента
        for (const auto &elem: elPool.elements) {
            elemUs[elem.ind][5] = 0.0;
            elemUs[elem.ind][6] = 0.0;
            std::vector<double> centroid = getElementCentroid2D(elem, nodePool);
            for (const auto &edgeInd: elem.edgeIndexes) {
                Edge edge = edgePool.edges[edgeInd];
                if (edge.neighbourInd1 == elem.ind) {
                    // у первого соседа в эдже заданы ноды в порядке полодительного обхода и нормаль тоже
                    const auto nodeInElemeInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                          edge.nodeInd1);
                    int node_before_ind =
                            nodeInElemeInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                    nodeInElemeInd - 1);
                    //        auto itNode1 = std::find(elementNodes.begin(), elementNodes.end(), node1);
                    Node node_before = nodePool.getNode(node_before_ind);
                    elemUs[elem.ind][5] += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                    elemUs[elem.ind][6] += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
                } else {
                    // а вот для второго нужно умножать на -1 и в обратном порядке
                    const auto nodeInElemeInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                          edge.nodeInd2);
                    int node_before_ind =
                            nodeInElemeInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                    nodeInElemeInd - 1);
                    //        auto itNode1 = std::find(elementNodes.begin(), elementNodes.end(), node2);
                    Node node_before = nodePool.getNode(node_before_ind);
                    elemUs[elem.ind][5] -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                    elemUs[elem.ind][6] -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
                }
            }
        }
        ++iterations;
        if(iterations > MAX_ITERATIONS){
            break;
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



void writeVTU(const std::string& filename, const World& geometryWorld, const std::vector<std::vector<double>>& elemUs) {
    const NodePool& np = geometryWorld.getNodePool();
    const ElementPool& ep = geometryWorld.getElementPool();

    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << np.nodeCount << "\" NumberOfCells=\"" << ep.elCount << "\">\n";

    // Write points (nodes)
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& node : np.nodes) {
        file << "          " << node.x << " " << node.y << " " << node.z << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // Write cells (elements)
    file << "      <Cells>\n";

    // Connectivity
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : ep.elements) {
        for (const auto& nodeIndex : element.nodeIndexes) {
            file << "          " << nodeIndex << " ";
        }
        file << "\n";
    }
    file << "        </DataArray>\n";

    // Offsets
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (const auto& element : ep.elements) {
        offset += element.dim;
        file << "          " << offset << "\n";
    }
    file << "        </DataArray>\n";

    // Types
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (const auto& element : ep.elements) {
        if (element.dim == 3) {
            file << "          5\n"; // Triangle
        } else if (element.dim == 4) {
            file << "          9\n"; // Quadrilateral
        }
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // Write solution data (elemUs)
    file << "      <CellData Scalars=\"elemUs\">\n";
    file << R"(        <DataArray type="Float64" NumberOfComponents=")" << elemUs[0].size() << "\" Name=\"elemUs\" format=\"ascii\">\n";
    for (const auto& U : elemUs) {
        for (const auto& value : U) {
            file << "          " << value << " ";
        }
        file << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </CellData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
}


void writeUnstructuredVTU(const std::string& filename,
                          const std::vector<std::array<double, 3>>& points,
                          const std::vector<std::vector<int>>& connectivity,
                          const std::vector<int>& cellTypes,
                          const std::vector<double>& cellValues) {
    std::ofstream file(filename);
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "<UnstructuredGrid>\n";
    file << "  <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << connectivity.size() << "\">\n";

    // Write Points
    file << "    <Points>\n";
    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& p : points) {
        file << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    file << "      </DataArray>\n";
    file << "    </Points>\n";

    // Write Cells
    file << "    <Cells>\n";

    // Connectivity
    file << "      <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& cell : connectivity) {
        for (int index : cell) {
            file << index << " ";
        }
        file << "\n";
    }
    file << "      </DataArray>\n";

    // Offsets
    file << "      <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (const auto& cell : connectivity) {
        offset += cell.size();
        file << offset << "\n";
    }
    file << "      </DataArray>\n";

    // Cell Types
    file << "      <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int type : cellTypes) {
        file << type << "\n";
    }
    file << "      </DataArray>\n";

    file << "    </Cells>\n";

    // Write Cell Data (Scalar Values)
    file << "    <CellData Scalars=\"Values\">\n";
    file << "      <DataArray type=\"Float64\" Name=\"Value\" format=\"ascii\">\n";
    for (double value : cellValues) {
        file << value << "\n";
    }
    file << "      </DataArray>\n";
    file << "    </CellData>\n";

    file << "  </Piece>\n";
    file << "</UnstructuredGrid>\n";
    file << "</VTKFile>\n";
}


