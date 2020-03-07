//
//  BnC.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 28/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ck_route/RouteMM.hpp"
#include "../../include/ck_route/CutBase.hpp"


rmm::BnC::BnC(rut::Problem* prob, std::string logPath, std::vector<std::string> cut_names, bool isTightenModel, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd, TimeTracker* tt) : RouteMM(prob, logPath, 'I', isTightenModel) {
    std::vector<CutBase*> cuts = get_cutInstances(cut_names, cutManagerType, isLocalCutAdd);
    cc = new CutComposer(prob, cuts, env, x_ij, logPath, tt);
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

int rmm::BnC::getNumGenCuts() {
    return cc->numGenCuts;
}

double rmm::BnC::getTime4Sep() {
    return cc->time4Sep;
}

double rmm::BnC::getTime4FDV() {
    return cc->time4FDV;
}

int rmm::BnC::getNum4Sep() {
    return cc->num4Sep;
}

void rmm::BnC::add_detectedCuts2MM() {
    for (CutBase *c: cc->cuts) {
        cplexModel->add(c->get_detectedCuts(cc));
    }
}
