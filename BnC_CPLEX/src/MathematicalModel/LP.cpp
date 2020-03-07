//
//  LP.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 28/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ck_route/RouteMM.hpp"

rmm::LP::LP(rut::Problem *prob, std::string logPath, bool isTightenModel) : RouteMM(prob, logPath, 'L', isTightenModel) {
    logStream = new std::ofstream (logPath.c_str(), std::ios_base::app);
    cplex->setOut(*logStream);
}

rmm::LP::~LP() {
    delete logStream;
}
