//
// Created by Иван on 10/21/2024.
//

#include "GasSolver2D.h"

void solverHLL2D(const World& world){
    NodePool np = world.getNodePool();
    ElementPool ep = world.getElementPool();
    std::vector<std::vector<double>> states(np.nodeCount, std::vector<double>(5, 0.0));
    /*
     *     states[i] <---> np.nodes[i]
     *
     * */

    //test1:   [0, 1] x [0, 1]
    /*   quadrant            rho        u         v         p
     *   x>0.5, y>0.5:          0.5313     0.0       0.0       0.4
     *   x<0.5, y>0.5:          1.0        0.7276    0.0       1.0
     *   x<0.5, y<0.5:          0.8        0.0       0.0       1.0
     *   x>0.5, y<0.5:          1.0        0.0       0.7276    1.0
     *   FREE-FLOW boundary conditions
     * */
    //printf("elCount = %d", ep.elCount);

    //initializer
    double gam_hcr = 3.0/2.0;
    double w = 0.0;
    double rho1 = 0.5313; double rho2 = 1.0; double rho3 = 0.8; double rho4 = 1.0;
    double u1 = 0.0; double u2 = 0.7276; double u3 = 0.0; double u4 = 0.0;
    double v1 = 0.0; double v2 = 0.0; double v3 = 0.0; double v4 = 0.7276;
    double p1 = 0.4; double p2 = 1.0; double p3 = 1.0; double p4 = 1.0;
    for(int i = 0; i < np.nodeCount; ++i){
        Node tmp_node = np.nodes[i];
        if(tmp_node.x > 0.5 && tmp_node.y > 0.5){
           // states[i] = state_from_primitive_vars(rho1, u1, v1, w, p1, gam_hcr);
        }
        else if(tmp_node.x < 0.5 && tmp_node.y > 0.5){
           // states[i] = state_from_primitive_vars(rho2, u2, v2, w, p2, gam_hcr);
        }
        else if(tmp_node.x < 0.5 && tmp_node.y < 0.5){
           // states[i] = state_from_primitive_vars(rho3, u3, v3, w, p3, gam_hcr);
        }
        else{
           // states[i] = state_from_primitive_vars(rho4, u4, v4, w, p4, gam_hcr);
        }
    }


    //solver
    for(int i = 0; i < ep.elCount; ++i){
        //пока делаем для n=4, для n=3 попозже...
        //std::vector<int> nodeInds = ep.elements[i].nodeIndexes;

    }

}