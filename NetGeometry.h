//
// Created by Иван on 10/21/2024.
//

#ifndef MAGNETTOPRJCT_NETGEOMETRY_H
#define MAGNETTOPRJCT_NETGEOMETRY_H

#include <vector>

//point
class Node{
public:
    double x;
    double y;
    double z;
};

class Element{
public:
    std::vector<double> nodes;
};


#endif //MAGNETTOPRJCT_NETGEOMETRY_H
