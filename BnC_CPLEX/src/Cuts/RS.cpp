//
//  RS.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 2/10/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ck_route/CutBase.hpp"


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

std::set<std::set<int>> solve_maxRS_constructCycle_RC(double **x_ij, rut::Problem *prob) {
    std::set<std::set<int>> outputSets;
    //
    int numNodes = (int) prob->N.size();
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
        }
        adjM_maxON.push_back(max_ON);
    }
    //
    bool visited[numNodes], isCycle;
    std::memset(visited, false, sizeof(visited));
    isCycle = false;
    //
    std::vector<int> path;
    
    int n0 = prob->o, n1;
    const int n2 = prob->d;
    path.push_back(n0);
    while (n0 != n2) {
        visited[n0] = true;
        n1 = adjM_maxON[n0];
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
    if (isCycle) {
        std::set<int> S1(path.begin(), path.end());
        outputSets.insert(S1);
    }
    return outputSets;
}

std::set<std::set<int>> solve_maxRS_constructCycle(CutComposer *cc, const IloCplex::Callback::Context &context) {
    std::set<std::set<int>> outputSets;
    rut::Problem *prob = cc->prob;
    //
    int numNodes = (int) prob->N.size();
    bool visited[numNodes], isCycle;
    std::memset(visited, false, sizeof(visited));
    isCycle = false;
    //
    std::vector<int> path;
    
    int n0 = prob->o, n1;
    const int n2 = prob->d;
    path.push_back(n0);
    while (n0 != n2) {
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
    if (isCycle) {
        std::set<int> S1(path.begin(), path.end());
        outputSets.insert(S1);
    }
    return outputSets;
}

RS_cut::RS_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd) : CutBase(ch_name, cutManagerType, isLocalCutAdd) {
    if (ch_name.substr(0, 1) == "e") {
        
    } else {
        assert(ch_name.substr(0, 1) == "h");
        separationAlgo = solve_maxRS_constructCycle;
    }
}

bool RS_cut::valid_subset(const std::set<int> &S1, rut::Problem *prob) {
    if ( S1.size() <= 1 ||
        generatedSets.find(S1) != generatedSets.end() ) {
        return false;
    }
    if (S1.find(prob->o) == S1.end()) {
        return false;
    }
    if (S1.find(prob->d) != S1.end()) {
        return false;
    }
    return true;
}

std::set<std::set<int>> RS_cut::solve_separationProb(CutComposer *cc, const IloCplex::Callback::Context &context) {
    rut::Problem *prob = cc->prob;
    std::set<std::set<int>> validSets;
    std::set<std::set<int>> outputSets = separationAlgo(cc, context);
    for (std::set<int> S1: outputSets) {
        if (valid_subset(S1, prob)) {
            double lhs = 0.0;
            for (int i: S1) {
                for (int j: S1) {
                    lhs += context.getRelaxationPoint(cc->x_ij[i][j]);
                }
            }
            int numFlowReturn = get_numFlowReturn(S1, prob);
            lhs -= ((int) S1.size() - (1 + numFlowReturn));
            if (lhs > 0) {
                validSets.insert(S1);
            }
        }
    }
    return validSets;
}

void RS_cut::set_LHS_Expr(IloExpr &lhs_expr, IloNumVar **x_ij, const std::set<int> &S1, rut::Problem *prob) {
    for (int i: S1) {
        for (int j: S1) {
            lhs_expr += x_ij[i][j];
        }
    }
    int numFlowReturn = get_numFlowReturn(S1, prob);
    lhs_expr -= ((int) S1.size() - (1 + numFlowReturn));
}

IloRangeArray RS_cut::get_cut_cnsts(double **x_ij, CutComposer *cc) {
    rut::Problem *prob = cc->prob;
    std::set<std::set<int>> validSets;
    std::set<std::set<int>> outputSets = solve_maxRS_constructCycle_RC(x_ij, prob);
    for (std::set<int> S1: outputSets) {
        if (valid_subset(S1, prob)) {
            double lhs = 0.0;
            for (int i: S1) {
                for (int j: S1) {
                    lhs += x_ij[i][j];
                }
            }
            int numFlowReturn = get_numFlowReturn(S1, prob);
            lhs -= ((int) S1.size() - (1 + numFlowReturn));
            if (lhs > 0) {
                validSets.insert(S1);
            }
        }
    }
    //
    char buf[2048];
    IloRangeArray cnsts(cc->env);
    IloExpr lhs_expr(cc->env);
    for (std::set<int> S1: validSets) {
        lhs_expr.clear();
        sprintf(buf, "%s(%d)", cut_name.c_str(), (int) generatedSets.size());
        generatedSets.insert(S1);
        set_LHS_Expr(lhs_expr, cc->x_ij, S1, cc->prob);
        cnsts.add(lhs_expr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    lhs_expr.end();
    return cnsts;
}

void RS_cut::add_cnsts2Model(const std::set<std::set<int>> &validSets,  CutComposer *cc, const IloCplex::Callback::Context &context) {
    for (std::set<int> S1: validSets) {
        generatedSets.insert(S1);
        IloExpr lhs_expr(cc->env);
        set_LHS_Expr(lhs_expr, cc->x_ij, S1, cc->prob);
        addUserCutwCust(context, lhs_expr);
        lhs_expr.end();
        cc->numGenCuts += 1;
    }
}

void RS_cut::add_cut(CutComposer *cc, const IloCplex::Callback::Context &context) {
    std::set<std::set<int>> validSets = solve_separationProb(cc, context);
    add_cnsts2Model(validSets, cc, context);
}

std::string RS_cut::add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context) {
    std::set<std::set<int>> validSets = solve_separationProb(cc, context);
    add_cnsts2Model(validSets, cc, context);
    std::string addedCuts;
    for (std::set<int> S1: validSets) {
        int numFlowReturn = get_numFlowReturn(S1, cc->prob);
        addedCuts += "(";
        for (int i: S1) {
            addedCuts += std::to_string(i) + "-";
        }
        addedCuts += ")" + std::to_string(numFlowReturn)+ ";";
    }
    //
    return addedCuts;
}

IloRangeArray RS_cut::get_detectedCuts(CutComposer *cc) {
    char buf[2048];
    IloRangeArray cnsts(cc->env);
    IloExpr lhs_expr(cc->env);
    for (std::set<int> S1: generatedSets) {
        if (addedSets2MM.find(S1) != addedSets2MM.end()) {
            continue;
        }
        addedSets2MM.insert(S1);
        lhs_expr.clear();
        sprintf(buf, "%s(%d)", cut_name.c_str(), (int) cnsts.getSize());
        set_LHS_Expr(lhs_expr, cc->x_ij, S1, cc->prob);
        cnsts.add(lhs_expr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    lhs_expr.end();
    return cnsts;
}
