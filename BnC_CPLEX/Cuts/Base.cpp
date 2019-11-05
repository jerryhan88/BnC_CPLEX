//
//  Base.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 13/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Base.hpp"

std::mutex mtx;

void CutBase::addUserCutwCust(const IloCplex::Callback::Context &context, IloExpr lhs_expr) {
    context.addUserCut(lhs_expr <= 0, cutManagerType, isLocalCutAdd);
}

CutComposer::CutComposer(Problem *prob, std::vector<CutBase*> &cuts, IloEnv &env, IloNumVar **x_ij, std::string logPath) {
    this->prob = prob;
    for (CutBase *c: cuts) {
        this->cuts.push_back(c);
    }
    this->env = env;
    this->x_ij = x_ij;
    this->logPath = logPath;
}

double** CutComposer::get_x_ij(const IloCplex::Callback::Context &context) {
    double **_x_ij = new double*[(*prob).N.size()];
    for (int i: (*prob).N) {
        _x_ij[i] = new double[(*prob).N.size()];
        for (int j: (*prob).N) {
            _x_ij[i][j] = context.getRelaxationPoint(x_ij[i][j]);
        }
    }
    return _x_ij;
}

void CutComposer::invoke (const IloCplex::Callback::Context &context) {
    if ( context.inRelaxation() ) {
        mtx.lock();
        double **_x_ij = get_x_ij(context);
        if (logPath == "") {
            for (CutBase *c: cuts) {
                c->add_cut(_x_ij, this, context);
            }
        } else {
            std::vector<std::string> violatedCnsts;
            for (CutBase *c: cuts) {
                std::string cnsts = c->add_cut_wLogging(_x_ij, this, context);
                violatedCnsts.push_back(cnsts);
            }
            //
            long nodecnt = context.getLongInfo(IloCplex::Callback::Context::Info::NodeCount);
            double objbst = context.getIncumbentObjective();
            std::string _objbst = objbst < 0.0 ? "-" : std::to_string(objbst);
            double objbnd = context.getRelaxationObjective();
            std::string _log(std::to_string(nodecnt) + ","
                             + _objbst + ","
                             + std::to_string(objbnd));
            for (std::string cnsts: violatedCnsts) {
                size_t nc = std::count(cnsts.begin(), cnsts.end(), ';');
                _log += "," + std::to_string(nc);
            }
            for (std::string cnsts: violatedCnsts) {
                _log += "," + cnsts;
            }
            char log[_log.size() + 1];
            std::strcpy(log, _log.c_str());
            appendRow(logPath, log);
        }
        mtx.unlock();
        for (int i: (*prob).N) {
            delete [] _x_ij[i];
        }
        delete [] _x_ij;
    }
}
