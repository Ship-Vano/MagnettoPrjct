//
// Created by Иван on 10/21/2024.
//

#include "NetGeometry.h"

#include <utility>

Node::Node(int xCoord, int yCoord, int zCoord) {
    x = xCoord;
    y = yCoord;
    z = zCoord;
}

ElementPool::ElementPool(int nodesPerElement, int elCnt, std::vector<Element> elements) {
    if(nodesPerElement == SQUARE_ELEMENT_NODE_COUNT){
        isSquare = true;
        isTriangular = false;
        for(int i = 0; i < elCnt; ++i){
            squareElements.emplace_back(elements[i]);
        }
    }
    else if(nodesPerElement == TRIANGULAR_ELEMENT_NODE_COUNT){
        isTriangular = true;
        isSquare = false;
        for(int i = 0; i < elCnt; ++i){
            triangularElements.emplace_back(elements[i]);
        }
    }
}

SquareElement::SquareElement(const Element& el) {
    nodeIndexes = el.nodeIndexes;
}

void World::setNodePool(const NodePool& np) {
    this->np = np;
}

NodePool World::getNodePool() {
    return np;
}

void World::setElementPool(const ElementPool &ep) {
    this->ep = ep;
}

ElementPool World::getElementPool() {
    return ep;
}

World::World(std::string filename) {

}


