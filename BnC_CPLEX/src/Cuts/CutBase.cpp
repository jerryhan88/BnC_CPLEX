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


int getNextNodeByFlow(int n0, CutComposer *cc, const IloCplex::Callback::Context &context) {
    int n1 = -1;
    double max_OEV = -1.0;
    for (int i = 0; i < cc->prob->N.size(); i++) {
        if (n0 == i) {
            continue;
        }
        if (!(cc->bool_x_ij[n0][i])) {
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

int get_numP(std::set<int> S1, rut::Problem *prob) {
    int numP = 0;
    for (int i: S1) {
        if (prob->P.find(i) != prob->P.end()) {
            numP += 1;
        }
    }
    return numP;
}

double get_vS(const std::set<int> &S1, rut::Problem *prob) {
    double vS = 0.0;
    for (int n0: S1) {
        vS += prob->v_i[n0];
    }
    return vS;
}

double get_wS(const std::set<int> &S1, rut::Problem *prob) {
    double wS = 0.0;
    for (int n0: S1) {
        wS += prob->w_i[n0];
    }
    return wS;
}

double get_LHS_CA(const std::set<int> &S, double **x_ij, rut::Problem *prob) {
    double lhs = 0.0;
    for (int i: S) {
        for (int j: S) {
            lhs += x_ij[i][j];
        }
    }
    int _vS, _wS, maxCapa;
    _vS = std::ceil(get_vS(S, prob) / prob->bv);
    _wS = std::ceil(get_wS(S, prob) / prob->bw);
    maxCapa = _vS > _wS ? _vS : _wS;
    lhs -= ((int) S.size() - maxCapa);
    return lhs;
}

int get_numFlowReturn(std::set<int> S1, rut::Problem *prob) {
    int numFlowReturn = 0;
    std::vector<int> included_newRR_indices;
    for (int i: S1) {
        if ( prob->R.find(i) != prob->R.end() &&
            i != prob->o ) {
            included_newRR_indices.push_back(i - prob->d);
        }
    }
    std::sort(included_newRR_indices.begin(), included_newRR_indices.end());
    for (int i = 0; i < included_newRR_indices.size(); i++) {
        if (i == 0) {
            numFlowReturn += included_newRR_indices[i] > 1 ? 1 : 0;
        } else {
            numFlowReturn += included_newRR_indices[i] - included_newRR_indices[i - 1] > 1 ? 1 : 0;
        }
    }
    return numFlowReturn;
}

bool get_cyclicPath(int ori, int dest, bool *visited, int numNodes, std::vector<int> &path,
                    CutComposer *cc, const IloCplex::Callback::Context &context) {
    std::memset(visited, false, numNodes);
    bool isCycle = false;
    path.clear();
    int n0 = ori, n1;
    //
    path.push_back(n0);
    while (n0 != dest) {
        visited[n0] = true;
        n1 = getNextNodeByFlow(n0, cc, context);
        if (n1 == -1) {
            break;
        }
        path.push_back(n1);
        if (visited[n1]) {
            isCycle = true;
            break;
        }
        n0 = n1;
    }
    return isCycle;
}

std::set<std::set<int>> get_maxCycles(const std::vector<int> &origins,
                   CutComposer *cc, const IloCplex::Callback::Context &context) {
    std::set<std::set<int>> outputSets;
    rut::Problem *prob = cc->prob;
    int numNodes = (int) prob->N.size();
    bool visited[numNodes], isCycle;
    const int n2 = prob->d;
    std::vector<int> path;
    for (int n0: origins) {
        isCycle = get_cyclicPath(n0, n2, visited, numNodes, path, cc, context);
        if (isCycle) {
            std::set<int> S1(path.begin(), path.end());
            outputSets.insert(S1);
        }
    }
    return outputSets;
}

std::set<std::set<edge>> get_infeasiblePaths(CutComposer *cc, const IloCplex::Callback::Context &context) {
    std::set<std::set<edge>> outputSets;
    rut::Problem *prob = cc->prob;
    //
    double *al_i = prob->al_i;
    double *be_i = prob->be_i;
    double **t_ij = prob->t_ij;
    int numNodes = (int) prob->N.size();
    //
    bool visited[numNodes], isCycle, lowLHS;
    bool is_TW_violated, is_DT_violated, is_RS_violated;
    double erest_deptTime, erest_arrvTime, tt, LHS;
    for (int i = 0; i < prob->_R.size() - 1; i++) {
        isCycle = false;
        lowLHS = false;
        std::memset(visited, false, sizeof(visited));
        is_TW_violated = false;
        is_DT_violated = false;
        is_RS_violated = false;
        //
        std::vector<edge> path;
        int n0 = prob->_R[i], n1;
        const int n2 = prob->_R[i + 1];
        erest_deptTime = al_i[n0];
        tt = 0.0;
        LHS = 1.0;
        while (n0 != n2) {
            visited[n0] = true;
            n1 = getNextNodeByFlow(n0, cc, context);
            if (n1 == -1) {
                break;
            }
            if (visited[n1]) {
                isCycle = true;
                break;
            }
            path.push_back(edge {n0, n1});
            LHS += (context.getRelaxationPoint(cc->x_ij[n0][n1]) - 1);
            if (LHS <= 0) {
                lowLHS = true;
                break;
            }
            if (n1 == prob->d && n2 != prob->d) {
                break;
            }
            tt += t_ij[n0][n1];
            erest_arrvTime = erest_deptTime + t_ij[n0][n1];
            if (be_i[n1] < erest_arrvTime) {
                is_TW_violated = true;
                break;
            }
            if (prob->bu < tt) {
                is_DT_violated = true;
                break;
            }
            if (prob->LR[i].find(n1) != prob->LR[i].end()) {
                is_RS_violated = true;
                break;
            }
            erest_deptTime = erest_arrvTime > al_i[n1] ? erest_arrvTime : al_i[n1];
            n0 = n1;
        }
        if (!isCycle && !lowLHS) {
            if (is_TW_violated || is_DT_violated || is_RS_violated) {
                std::set<edge> S1(path.begin(), path.end());
                outputSets.insert(S1);
            }
        }
    }
    //
    return outputSets;
}


void CutBase::addUserCutwCust(const IloCplex::Callback::Context &context, IloExpr lhs_expr) {
    context.addUserCut(lhs_expr <= 0, cutManagerType, isLocalCutAdd);
}

void CutComposer::invoke(const IloCplex::Callback::Context &context) {
    if ( context.inRelaxation() ) {
        cut_mtx.lock();
        //
        double relaxedVal = context.getRelaxationObjective();
        double objbst = context.getIncumbentObjective();
        if (bestRelVal < objbst) {
            bestRelVal = DBL_MAX;
        }
        if (relaxedVal < bestRelVal) {
            bestRelVal = relaxedVal;
            if (logPath == "") {
                for (CutBase *c: cuts) {
                    c->add_cut(this, context);
                }
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
