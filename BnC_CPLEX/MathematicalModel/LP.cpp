//
//  LP.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 28/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "BaseMM.hpp"

LP::LP(Problem *prob, std::string logPath, bool isTightenModel) : BaseMM(prob, logPath, 'L', isTightenModel) {
    logStream = new std::ofstream (logPath.c_str(), std::ios_base::app);
    cplex->setOut(*logStream);
}

LP::~LP() {
    delete logStream;
}
