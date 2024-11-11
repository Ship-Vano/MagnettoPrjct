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
#include <cmath>

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
    std::vector<double> normalVector; // компоненты вектора нормали
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

    Node getNode(int ind);
};

// Набор рёбер
class EdgePool{
public:
    int edgeCount;
    std::vector<Edge> edges;
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

double areaCalc(const Element& poly, const NodePool& nPool);

#endif // MAGNETTOPRJCT_NETGEOMETRY_H
