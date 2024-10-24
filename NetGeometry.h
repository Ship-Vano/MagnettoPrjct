//
// Created by Иван on 10/21/2024.
//

#ifndef MAGNETTOPRJCT_NETGEOMETRY_H
#define MAGNETTOPRJCT_NETGEOMETRY_H

#include "NamesNConstants.h"
#include <vector>
#include <string>
#include <fstream>

// Point
class Node{
public:
    double x;
    double y;
    double z;
    Node(int xCoord, int yCoord, int zCoord);
};

// Mesh elements
class Element{
public:
    std::vector<double> nodeIndexes;
};

class SquareElement: public Element{
public:
 const int nodeSize = 4;

 explicit SquareElement(const Element& el);
};

class TriangularElement: public Element{
public:
    const int nodeSize = 3;
};

// Node pool
class NodePool{
public:
    int nodeCount;
    std::vector<Node> nodes;
};

// Element pool
class ElementPool{
public:
    int elCount;
    bool isSquare;
    bool isTriangular;
    std::vector<SquareElement> squareElements;
    std::vector<TriangularElement> triangularElements;
    ElementPool(int nodesPerElement, int elCnt, std::vector<Element> elements);
};

// Problem's World
class World{
private:
    NodePool np;
    ElementPool ep;
public:
    void setNodePool(const NodePool& np);
    NodePool getNodePool();
    void setElementPool(const ElementPool& ep);
    ElementPool getElementPool();
    World(std::string filename);
};

#endif //MAGNETTOPRJCT_NETGEOMETRY_H
