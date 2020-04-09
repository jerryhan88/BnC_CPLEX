//
//  ILP.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 27/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ck_route/Router.hpp"


ILOMIPINFOCALLBACK1(BnBLogger,
                    std::string, logPath) {
    long nodecnt = getNnodes();
    double objbst = 0.0;
    if ( hasIncumbent() ) {
        objbst = getIncumbentObjValue();
    }
    double objbnd = getBestObjValue();
    std::string _log(std::to_string(nodecnt) + ","
                     + std::to_string(objbst) + ","
                     + std::to_string(objbnd));
    char log[_log.size() + 1];
    std::strcpy(log, _log.c_str());
    appendRow(logPath, log);
}

rmm::ILP::ILP(rut::Problem *prob, TimeTracker* tt, unsigned long time_limit_sec, std::string logPath, std::string lpPath, bool isTightenModel) : RouteMM(prob, tt, time_limit_sec, logPath, lpPath, 'I', isTightenModel) {
    cplex->setOut(env.getNullStream());
    cplex->use(BnBLogger(env, logPath));
    if (logPath != "") {
        std::string _header("nCnt,bestObj,bestBound,note");
        char header[_header.size() + 1];
        std::strcpy(header, _header.c_str());
        createCSV(this->logPath, header);
    }
}

