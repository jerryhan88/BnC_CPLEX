//
//  RC.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 28/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ck_route/Router.hpp"
#include "../../include/ck_route/CutBase.hpp"

rut::Solution* rmm::RC::solve() {
    double** _x_ij = new double *[(*prob).N.size()];
    double* _u_i = new double[(*prob).N.size()];
    for (int i: (*prob).N) {
        _x_ij[i] = new double[(*prob).N.size()];
    }
    int counter = 0;
    while (true) {
        cplex->solve();
        std::string _log(std::to_string(counter) + "," + std::to_string(cplex->getObjValue()));
        get_x_ij(_x_ij);
        get_u_i(_u_i);
        //
        int numGenCuts = 0;
        for (CutBase *c: cc->cuts) {
            IloRangeArray cnsts = c->get_cut_cnsts(_x_ij, cc);
            cplexModel->add(cnsts);
            numGenCuts += cnsts.getSize();
            _log += "," + std::to_string(cnsts.getSize());
        }
        char log[_log.size() + 1];
        std::strcpy(log, _log.c_str());
        appendRow(logPath, log);
        if (numGenCuts == 0) {
            break;
        }
        counter++;
    }
    
    for (int i = 0; i < (*prob).N.size(); i++) {
        delete [] _x_ij[i];
    }
    delete [] _x_ij;
    delete [] _u_i;
    
    
    return nullptr;
}
