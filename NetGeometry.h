//
// Created by Иван on 10/21/2024.
//

#ifndef MAGNETTOPRJCT_NETGEOMETRY_H
#define MAGNETTOPRJCT_NETGEOMETRY_H

#include "NamesNConstants.h"
#include "LinOp.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <utility>
#include <iostream>
#include <cmath>
#include <cassert>
#include <map>
#include <unordered_map>
#include <unordered_set>

// Узел
class Node {
public:
    int ind; // порядковый номер
    double x;
    double y;
    double z;

    Node(int index, double xCoord, double yCoord, double zCoord);
};

// Ребро
class Edge{
public:
    int ind; // порядковый номер
    int nodeInd1; // номер первого узла
    int nodeInd2; // номер второго узла
    int neighbourInd1; //номер первого соседнего элемента
    int neighbourInd2; //номер второго соседнего элемента
    int orientation; // ориентация (по идее будет = 1 или = -1)
    double length;
    std::vector<double> normalVector; // компоненты вектора нормали
    Edge(int index, int node1, int node2, int neighbor1, int neighbor2,
         int orient, double len, const std::vector<double>& normalVec);
};

// Элемент
class Element {
public:
    int ind; // порядковый номер
    int dim; // размерность (кол-во узлов)
    std::vector<int> nodeIndexes; //номера узлов
    std::vector<int> edgeIndexes; //номера рёбер
    double area = 0.; //площадь
    Element(const int index, const std::vector<int> &nIndexes, int size);
};

// Набор узлов
class NodePool {
public:
    int nodeCount;
    std::vector<Node> nodes;

    NodePool(int size, const std::vector<Node>& nodeVec);
    NodePool() : nodeCount(0), nodes() {} // Default constructor

    Node getNode(int ind) const;
};

// Набор элементов
class ElementPool {
public:
    int elCount;
    bool isSquare;
    bool isTriangular;
    std::vector<Element> elements;

    ElementPool(int nodesPerElement, int elCnt, const std::vector<Element>& elements);
    ElementPool() : elCount(0), isSquare(false), isTriangular(false), elements() {} // Default constructor
};

// Набор рёбер
class EdgePool{
public:
    int edgeCount;
    std::vector<Edge> edges;
    EdgePool(int size, const std::vector<Edge>& edgeVec);
    EdgePool(const NodePool& np, ElementPool& ep);
    EdgePool() : edgeCount(0), edges() {} // Default constructor
};

class NeighbourService {
private:
    std::unordered_map<int, std::unordered_set<int>> nodeToElements; // Nodes -> Elements
    std::unordered_map<int, std::vector<int>> nodeToEdgesMap;       //Nodes -> Edges
    std::unordered_map<int, std::unordered_set<int>> edgeToElements; // Edges -> Elements
    std::unordered_map<int, std::unordered_set<int>> elementToElements; // Elements -> Elements
    std::unordered_map<int, std::unordered_set<int>> elementToEdges; // Elements -> Edges
    std::unordered_map<int, std::unordered_set<int>> edgeToElementsMap; // Edges -> Elements

public:
    NeighbourService(const NodePool& np, const ElementPool& ep, const EdgePool& edgePool);
    std::vector<int> getEdgeNeighborsOfNode(int nodeIndex) const;
    std::unordered_set<int> getNodeNeighbours(int nodeIndex) const;
    std::unordered_set<int> getEdgeNeighbours(int edgeIndex) const;
    std::unordered_set<int> getElementNeighbours(int elementIndex) const;
    std::unordered_set<int> getEdgesOfElement(int elementIndex) const;
    std::unordered_set<int> getElementsOfEdge(int edgeIndex) const;
    int findEdgeByNodes(int node1Index, int node2Index, const EdgePool& edgePool) const;

    void displayNeighbours() const;
};


// Problem's World
class World {
private:
    NodePool np;
    ElementPool ep;
    EdgePool edgp;
    NeighbourService ns;
public:
    void setNodePool(const NodePool& np);
    NodePool getNodePool() const;
    void setElementPool(const ElementPool& ep);
    ElementPool getElementPool() const;
    EdgePool getEdgePool() const;
    void setEdgePool(const EdgePool& edgp);
    World(const std::string& fileName);
    void display() const;
    NeighbourService& getNeighbourService();
};

double areaCalc(const Element& poly, const NodePool& nPool);

std::vector<double> getElementCentroid2D(const Element& poly, const NodePool& nPool);

std::vector<double> getMidPoint2D(const int nodeInd1, const int nodeInd2, const NodePool& nPool);

double getDistance(const int nodeInd1, const int nodeInd2, const NodePool& nPool);

void setNeighbourEdge(Element& el, const int edgeInd);
#endif // MAGNETTOPRJCT_NETGEOMETRY_H
