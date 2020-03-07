//
//  CutBase.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 1/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ck_route/CutBase.hpp"

#include <random>

std::mutex cut_mtx;

#define REL_VALS_LIMIT 32

int getNextNodeByFlow(int n0, CutComposer *cc, const IloCplex::Callback::Context &context) {
    int n1 = -1;
    double max_OEV = -1.0;
    for (int i = 0; i < cc->prob->N.size(); i++) {
        if (n0 == i) {
            continue;
        }
        double oev = context.getRelaxationPoint(cc->x_ij[n0][i]);
        if (oev > 0.0 && max_OEV < oev) {
            n1 = i;
            max_OEV = oev;
        }
        if (oev > 0.5) {
            break;
        }
    }
    return n1;
}

void CutBase::addUserCutwCust(const IloCplex::Callback::Context &context, IloExpr lhs_expr) {
    context.addUserCut(lhs_expr <= 0, cutManagerType, isLocalCutAdd);
}

CutComposer::CutComposer(rut::Problem *prob, std::vector<CutBase*> &cuts, IloEnv &env, IloNumVar **x_ij, std::string logPath, TimeTracker* tt) {
    this->prob = prob;
    for (CutBase *c: cuts) {
        this->cuts.push_back(c);
    }
    this->env = env;
    this->x_ij = x_ij;
    this->logPath = logPath;
    this->tt = tt;
    //
    _x_ij = new IloNum*[(*prob).N.size()];
    for (int i: (*prob).N) {
        _x_ij[i] = new IloNum[(*prob).N.size()];
    }
}

void CutComposer::invoke (const IloCplex::Callback::Context &context) {
    if ( context.inRelaxation() ) {
        cut_mtx.lock();
        //
        double relaxedVal = context.getRelaxationObjective();
        double objbst = context.getIncumbentObjective();
        if (bestRelVal < objbst) {
            bestRelVal = DBL_MAX;
        }
//        if (relaxedVals.find(relaxedVal) == relaxedVals.end()) {
//            if (relaxedVals.size() == REL_VALS_LIMIT) {
//                relaxedVals.clear();
//            }
        if (relaxedVal < bestRelVal) {
            bestRelVal = relaxedVal;

            
//            double timeRecored = tt->get_elipsedTimeCPU();
//            time4FDV += tt->get_elipsedTimeCPU() - timeRecored;
//            timeRecored = tt->get_elipsedTimeCPU();
            if (logPath == "") {
                for (CutBase *c: cuts) {
                    c->add_cut(this, context);
                }
//                time4Sep += tt->get_elipsedTimeCPU() - timeRecored;
                num4Sep += 1;
            } else {
                double timeRecored = tt->get_elapsedTimeCPU();
                std::vector<std::string> violatedCnsts;
                for (CutBase *c: cuts) {
                    std::string cnsts = c->add_cut_wLogging(this, context);
                    violatedCnsts.push_back(cnsts);
                }
                time4Sep += tt->get_elapsedTimeCPU() - timeRecored;
                num4Sep += 1;
                //
                long nodecnt = context.getLongInfo(IloCplex::Callback::Context::Info::NodeCount);
                
                std::string _objbst = objbst < 0.0 ? "-" : std::to_string(objbst);
                std::string _log(std::to_string(nodecnt) + ","
                                 + _objbst + ","
                                 + std::to_string(relaxedVal));
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
        }
        //
        cut_mtx.unlock();
    }
}
