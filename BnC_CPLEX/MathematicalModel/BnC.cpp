//
//  BnC.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 28/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "BaseMM.hpp"
#include "../Cuts/Base.hpp"


BnC::BnC(Problem* prob, std::string logPath, std::vector<std::string> cut_names, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd) : BaseMM(prob, logPath, 'I', true) {
    std::vector<CutBase*> cuts = get_cutInstances(cut_names, cutManagerType, isLocalCutAdd);
    cc = new CutComposer(prob, cuts, env, x_ij, logPath);
    CPXLONG contextMask = 0;
    contextMask |= IloCplex::Callback::Context::Id::Relaxation;
    cplex->use(cc, contextMask);
    //
    cplex->setOut(env.getNullStream());
    if (logPath != "") {
        std::string _header("nCnt,bestObj,bestBound");
        for (CutBase *c: cc->cuts) {
            _header += "," + c->cut_name;
        }
        _header += ",note";
        char header[_header.size() + 1];
        std::strcpy(header, _header.c_str());
        createCSV(logPath, header);
    }
}

int BnC::getNumGenCuts() {
    return cc->numGenCuts;
}
