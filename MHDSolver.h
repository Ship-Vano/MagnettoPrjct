//
// Created by Иван on 4/29/2024.
//

#ifndef MAGNETTOPRJCT_MHDSOLVER_H
#define MAGNETTOPRJCT_MHDSOLVER_H

#include "MHDProblem.h"
#include "LinOp.h"
#include "FileIO.h"

std::vector<double> state_from_primitive_vars(const double &rho, const double &u, const double &v, const double &w, const double &p, const double &Bx, const double &By, const double &Bz, const double &gam_hcr);

bool HLLScheme(const MHDProblem &problem, const std::string &filename ="HLLScheme");

bool HLLCScheme(const MHDProblem &problem, const std::string &filename ="HLLCScheme");

bool HLLDScheme(const MHDProblem &problem, const std::string &filename ="HLLDScheme");

#endif //MAGNETTOPRJCT_MHDSOLVER_H
