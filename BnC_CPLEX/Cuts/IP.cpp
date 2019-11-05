//
//  IP.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 15/10/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Base.hpp"




std::set<std::set<edge>> solve_allNodes(double** x_ij, Problem *prob) {
    std::set<std::set<edge>> outputSets;
    //
    double* al_i = (*prob).al_i;
    double* be_i = (*prob).be_i;
    double** t_ij = (*prob).t_ij;
    int numNodes = (int) (*prob).N.size();
    std::vector<int> adjM_maxON;
    for (int i = 0; i < numNodes; i++) {
        int max_ON = -1;
        double max_OEV = -1.0;
        for (int j = 0; j < numNodes; j++) {
            if (i == j) {
                continue;
            }
            if (x_ij[i][j] > 0 && max_OEV < x_ij[i][j]){
                max_ON = j;
                max_OEV = x_ij[i][j];
            }
//            std::cout << "x[" << i << "][" << j << "]:" << x_ij[i][j] << std::endl;
        }
        adjM_maxON.push_back(max_ON);
    }
    //
    bool visited[numNodes], isCycle;
    bool is_TW_violated, is_DT_violated, is_RS_violated;
    double erest_deptTime, erest_arrvTime, tt;
    for (int i = 0; i < (*prob)._S.size() - 1; i++) {
        isCycle = false;
        std::memset(visited, false, sizeof(visited));
        is_TW_violated = false;
        is_DT_violated = false;
        is_RS_violated = false;
        //
        std::vector<edge> path;
        int n0 = (*prob)._S[i], n1;
        const int n2 = (*prob)._S[i + 1];
        erest_deptTime = al_i[n0];
        tt = 0.0;
        while (n0 != n2) {
            visited[n0] = true;
            n1 = adjM_maxON[n0];
            if (n1 == -1) {
                break;
            }
            if (visited[n1]) {
                isCycle = true;
                break;
            }
            path.push_back(edge {n0, n1});
            if (n1 == (*prob).d && n2 != (*prob).d) {
                break;
            }
            tt += t_ij[n0][n1];
            erest_arrvTime = erest_deptTime + t_ij[n0][n1];
            if (be_i[n1] < erest_arrvTime) {
                is_TW_violated = true;
                break;
            }
            if ((*prob).bu < tt) {
                is_DT_violated = true;
                break;
            }
            if ((*prob).LS[i].find(n1) != (*prob).LS[i].end()) {
                is_RS_violated = true;
                break;
            }
            erest_deptTime = erest_arrvTime > al_i[n1] ? erest_arrvTime : al_i[n1];
            n0 = n1;
        }
        if (!isCycle) {
            if (is_TW_violated || is_DT_violated || is_RS_violated) {
                std::set<edge> S1(path.begin(), path.end());
                outputSets.insert(S1);
            }
        }
    }
    //
    for (int k: (*prob).K) {
        isCycle = false;
        std::memset(visited, false, sizeof(visited));
        is_TW_violated = false;
        is_DT_violated = false;
        
        std::vector<edge> path;
        int n0 = (*prob).h_k[k], n1;
        const int n2 = (*prob).n_k[k];
        erest_deptTime = al_i[n0];
        tt = 0.0;
        while (n0 != n2) {
            visited[n0] = true;
            n1 = adjM_maxON[n0];
            if (n1 == -1) {
                break;
            }
            if (visited[n1]) {
                isCycle = true;
                break;
            }
            path.push_back(edge {n0, n1});
            if (n1 == (*prob).d) {
                break;
            }
            tt += t_ij[n0][n1];
            erest_arrvTime = erest_deptTime + t_ij[n0][n1];
            if (be_i[n1] < erest_arrvTime) {
                is_TW_violated = true;
                break;
            }
            if ((*prob).bu < tt) {
                is_DT_violated = true;
                break;
            }
            erest_deptTime = erest_arrvTime > al_i[n1] ? erest_arrvTime : al_i[n1];
            n0 = n1;
        }
        if (!isCycle) {
            if (is_TW_violated || is_DT_violated) {
                std::set<edge> S1(path.begin(), path.end());
                outputSets.insert(S1);
            }
        }
    }
    //
    return outputSets;
}

IP_cut::IP_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd) : CutBase(ch_name, cutManagerType, isLocalCutAdd) {
    separationAlgo = solve_allNodes;
}

bool IP_cut::valid_subset(const std::set<edge>& S1) {
    if (generatedSets.find(S1) != generatedSets.end() ) {
        return false;
    } else {
        return true;
    }
}

std::set< std::set<edge> > IP_cut::solve_separationProb(double** x_ij, CutComposer* cc) {
    Problem *prob = cc->prob;
    std::set<std::set<edge>> validSets;
    std::set<std::set<edge>> outputSets = separationAlgo(x_ij, prob);
    //
    for (std::set<edge> S1: outputSets) {
        if (valid_subset(S1)) {
            double lhs = 0.0;
            for (edge e: S1) {
                lhs += x_ij[e.first][e.second];
            }
            lhs -= ((int) S1.size() - 1);
            if (lhs > 0) {
                validSets.insert(S1);
            }
        }
    }
    return validSets;
}

void IP_cut::set_LHS_Expr(IloExpr& lhs_expr, IloNumVar** x_ij, const std::set<edge>& S1) {
    for (edge e: S1) {
        lhs_expr += x_ij[e.first][e.second];
    }
    lhs_expr -= ((int) S1.size() - 1);
}

IloRangeArray IP_cut::get_cut_cnsts(double** x_ij, CutComposer* cc) {
    std::set<std::set<edge>> validSets = solve_separationProb(x_ij, cc);
    char buf[2048];
    IloRangeArray cnsts(cc->env);
    for (std::set<edge> S1: validSets) {
        sprintf(buf, "%s(%d)", cut_name.c_str(), (int) generatedSets.size());
        generatedSets.insert(S1);
        IloExpr lhs_expr(cc->env);
        set_LHS_Expr(lhs_expr, cc->x_ij, S1);
        cnsts.add(lhs_expr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
        lhs_expr.end();
    }
    return cnsts;
}

void IP_cut::add_cnsts2Model(const std::set<std::set<edge>>& validSets, CutComposer* cc, const IloCplex::Callback::Context& context) {
    for (std::set<edge> S1: validSets) {
        generatedSets.insert(S1);
        IloExpr lhs_expr(cc->env);
        set_LHS_Expr(lhs_expr, cc->x_ij, S1);
        addUserCutwCust(context, lhs_expr);
        lhs_expr.end();
        cc->numGenCuts += 1;
    }
}

void IP_cut::add_cut(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context) {
    std::set<std::set<edge>> validSets = solve_separationProb(x_ij, cc);
    add_cnsts2Model(validSets, cc, context);
}

std::string IP_cut::add_cut_wLogging(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context) {
    std::set<std::set<edge>> validSets = solve_separationProb(x_ij, cc);
    add_cnsts2Model(validSets, cc, context);
    std::string addedCuts;
    for (std::set<edge> S1: validSets) {
        std::map<int, int> _route;
        std::map<int, int> _routeRev;
        for (edge e: S1) {
            _route[e.first] = e.second;
            _routeRev[e.second] = e.first;
//            addedCuts += std::to_string(e.first) + "-" + std::to_string(e.second) + "#";
        }
        int n1 = (*S1.begin()).second;
        int n0;
        while (true) {
            if (_routeRev.find(n1) != _routeRev.end()) {
                n1 = _routeRev[n1];
            } else {
                n0 = n1;
                break;
            }
        }
        addedCuts += "(" + std::to_string(n0);
        while (true) {
            if (_route.find(n0) != _route.end()) {
                n0 = _route[n0];
                addedCuts += "-" + std::to_string(n0);
            } else {
                break;
            }
        }
        addedCuts += ");";
    }
    return addedCuts;
}
