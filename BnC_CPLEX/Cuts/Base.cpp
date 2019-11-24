//
//  Base.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 13/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Base.hpp"
#include <random>


std::mutex mtx;

void CutBase::addUserCutwCust(const IloCplex::Callback::Context &context, IloExpr lhs_expr) {
    context.addUserCut(lhs_expr <= 0, cutManagerType, isLocalCutAdd);
}

CutComposer::CutComposer(Problem *prob, std::vector<CutBase*> &cuts, IloEnv &env, IloNumVarArray* x_ij, std::string logPath, TimeTracker* tt) {
    this->prob = prob;
    for (CutBase *c: cuts) {
        this->cuts.push_back(c);
    }
    this->env = env;
    this->x_ij = x_ij;
    this->logPath = logPath;
    this->tt = tt;
    //
    _x_ij = IloArray<IloNumArray>(this->env, (*prob).N.size());
//    _x_ij = new IloNum*[(*prob).N.size()];
    for (int i: (*prob).N) {
        _x_ij[i] = IloNumArray(this->env, (*prob).N.size());
        
//        _x_ij[i] = new IloNum[(*prob).N.size()];
    }
}

void CutComposer::set_x_ij(const IloCplex::Callback::Context &context) {
    for (int i: (*prob).N) {
        
        context.getRelaxationPoint(x_ij[i], _x_ij[i]);
        
        
//        for (int j: (*prob).N) {
//            _x_ij[i][j] = context.getRelaxationPoint(x_ij[i][j]);
//        }
    }
}

void CutComposer::invoke (const IloCplex::Callback::Context &context) {
    if ( context.inRelaxation() ) {
        mtx.lock();
//        double **_x_ij = get_x_ij(context);
        double timeRecored = tt->get_elipsedTimeCPU();
        set_x_ij(context);
        time4FDV += tt->get_elipsedTimeCPU() - timeRecored;
        timeRecored = tt->get_elipsedTimeCPU();
        if (logPath == "") {
            for (CutBase *c: cuts) {
                c->add_cut(_x_ij, this, context);
            }
            time4Sep += tt->get_elipsedTimeCPU() - timeRecored;
            num4Sep += 1;
        } else {
            std::vector<std::string> violatedCnsts;
            for (CutBase *c: cuts) {
                std::string cnsts = c->add_cut_wLogging(_x_ij, this, context);
                violatedCnsts.push_back(cnsts);
            }
            time4Sep += tt->get_elipsedTimeCPU() - timeRecored;
            num4Sep += 1;
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
//        for (int i: (*prob).N) {
//            delete [] _x_ij[i];
//        }
//        delete [] _x_ij;
        mtx.unlock();
    }
}
