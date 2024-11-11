//
// Created by Иван on 10/21/2024.
//

#include <cassert>
#include "NetGeometry.h"

Node::Node(int index, double xCoord, double yCoord, double zCoord)
        : ind(index), x(xCoord), y(yCoord), z(zCoord) {}

Element::Element(const int index, const std::vector<int> &nIndexes, int size)
        : ind(index), nodeIndexes(nIndexes), dim(size) {
}

double areaCalc(const Element& poly, const NodePool& nPool) {
    int dim = poly.dim;
    std::vector<Node> polyNodes;
    for(int i = 0; i < dim; ++i){
        polyNodes.push_back(nPool.nodes[poly.nodeIndexes[i]]);
    }
    double res = 0.0;
    double iSum = 0.0;
    double jSum = 0.0;
    double kSum = 0.0;
    for(int k = 1; k < dim-1; ++k){
        double yk_y1 = polyNodes[k].y - polyNodes[0].y;  //y_k - y_1
        double zk1_z1 = polyNodes[k+1].z - polyNodes[0].z; // z_{k+1} - z_1
        double zk_z1 = polyNodes[k].z - polyNodes[0].z;
        double yk1_y1 = polyNodes[k+1].y - polyNodes[0].y;
        double xk1_x1 = polyNodes[k+1].x - polyNodes[0].x;
        double xk_x1 = polyNodes[k].x - polyNodes[0].x;
        iSum += (yk_y1 * zk1_z1 - zk_z1 * yk1_y1);
        jSum += (zk_z1 * xk1_x1 - xk_x1 * zk1_z1);
        kSum += (xk_x1 * yk1_y1 - yk_y1 * xk1_x1);
    }
    res = 0.5 * std::sqrt(iSum*iSum + jSum*jSum + kSum*kSum);
    return res;
}

NodePool::NodePool(int size, const std::vector<Node>& nodeVec)
        : nodeCount(size), nodes(nodeVec) {}

Node NodePool::getNode(int ind) {
    return nodes[ind];
}

ElementPool::ElementPool(int nodesPerElement, int elCnt, const std::vector<Element>& elems)
        : elCount(elCnt), isSquare(nodesPerElement == SQUARE_ELEMENT_NODE_COUNT),
          isTriangular(nodesPerElement == TRIANGULAR_ELEMENT_NODE_COUNT), elements(elems) {}

void World::setNodePool(const NodePool& np) {
    this->np = np;
}

NodePool World::getNodePool() const {
    return np;
}

void World::setElementPool(const ElementPool& ep) {
    this->ep = ep;
}

ElementPool World::getElementPool() const {
    return ep;
}

World::World(const std::string& fileName) : np(), ep() {
    std::ifstream file(fileName);
    assert(file.is_open());

    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::string tmp_line;

    while (std::getline(file, tmp_line)) {
        if (tmp_line == "$Nodes") {
            std::getline(file, tmp_line);  // Size
            std::getline(file, tmp_line);  // First enter
            while (tmp_line != "$EndNodes") {
                int ind;
                double x, y, z;
                std::istringstream ss(tmp_line);
                ss >> ind >> x >> y >> z;  // Read all values in one go
                nodes.emplace_back(ind, x, y, z);
                std::getline(file, tmp_line);
            }
        } else if (tmp_line == "$Elements") {
            std::getline(file, tmp_line); // Size
            std::getline(file, tmp_line); // First enter
            while (tmp_line != "$EndElements") {
                int ind, count;
                std::istringstream ss(tmp_line);
                ss >> ind >> count; // Read index and count

                std::vector<int> indexes(count);
                for (int i = 0; i < count; ++i) {
                    ss >> indexes[i]; // Directly read into vector
                }
                elements.emplace_back(ind, indexes, count);
                std::getline(file, tmp_line);
            }
        }
    }

    np = NodePool(nodes.size(), nodes);
    ep = ElementPool(nodes[0].ind, elements.size(), elements); // Assuming nodes[0].ind is the nodesPerElement

    for(int i = 0; i < ep.elCount; ++i){
        ep.elements[i].area = areaCalc(ep.elements[i], np);
    }

}

void World::display() const {
    std::cout << "Node Pool:" << std::endl;
    std::cout << "Total Nodes: " << np.nodeCount << std::endl;
    for (const auto& node : np.nodes) {
        std::cout << "Node Index: " << node.ind << ", Coordinates: ("
                  << node.x << ", " << node.y << ", " << node.z << ")" << std::endl;
    }

    std::cout << "Element Pool:" << std::endl;
    std::cout << "Total Elements: " << ep.elCount << std::endl;
    for (const auto& element : ep.elements) {
        std::cout << "Element with Node Size: " << element.dim << ", Area = " << element.area<<", Node Indices: ";
        for (const auto& index : element.nodeIndexes) {
            std::cout << index << " ";
        }
        std::cout << std::endl;
    }
}