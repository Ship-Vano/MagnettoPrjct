//
// Created by Иван on 10/21/2024.
//

#ifndef MAGNETTOPRJCT_NETGEOMETRY_H
#define MAGNETTOPRJCT_NETGEOMETRY_H

#include "NamesNConstants.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <utility>
#include <iostream>

// Point
class Node {
public:
    int ind;
    double x;
    double y;
    double z;

    Node(int index, double xCoord, double yCoord, double zCoord);
};

// Mesh elements
class Element {
public:
    int nodeSize; // Changed to non-const
    std::vector<int> nodeIndexes;

    Element(const std::vector<int>& indexes, int nodeSize);
};

// Node pool
class NodePool {
public:
    int nodeCount;
    std::vector<Node> nodes;

    NodePool(int size, const std::vector<Node>& nodeVec);
    NodePool() : nodeCount(0), nodes() {} // Default constructor
};

// Element pool
class ElementPool {
public:
    int elCount;
    bool isSquare;
    bool isTriangular;
    std::vector<Element> elements;

    ElementPool(int nodesPerElement, int elCnt, const std::vector<Element>& elements);
    ElementPool() : elCount(0), isSquare(false), isTriangular(false), elements() {} // Default constructor
};

// Problem's World
class World {
private:
    NodePool np;
    ElementPool ep;
public:
    void setNodePool(const NodePool& np);
    NodePool getNodePool() const;
    void setElementPool(const ElementPool& ep);
    ElementPool getElementPool() const;
    World(const std::string& fileName);
    void display() const;
};

#endif // MAGNETTOPRJCT_NETGEOMETRY_H
