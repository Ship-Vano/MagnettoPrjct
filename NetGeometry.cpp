//
// Created by Иван on 10/21/2024.
//

#include <cassert>
#include "NetGeometry.h"

Node::Node(int index, double xCoord, double yCoord, double zCoord)
        : ind(index), x(xCoord), y(yCoord), z(zCoord) {}

Element::Element(const std::vector<int>& indexes, int nodeSize)
        : nodeIndexes(indexes), nodeSize(nodeSize) {}

NodePool::NodePool(int size, const std::vector<Node>& nodeVec)
        : nodeCount(size), nodes(nodeVec) {}

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
                elements.emplace_back(indexes, count);
                std::getline(file, tmp_line);
            }
        }
    }

    np = NodePool(nodes.size(), nodes);
    ep = ElementPool(nodes[0].ind, elements.size(), elements); // Assuming nodes[0].ind is the nodesPerElement
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
        std::cout << "Element with Node Size: " << element.nodeSize << ", Node Indices: ";
        for (const auto& index : element.nodeIndexes) {
            std::cout << index << " ";
        }
        std::cout << std::endl;
    }
}