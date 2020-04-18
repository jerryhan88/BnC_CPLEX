//
//  SE.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 14/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//


#include "../../include/ck_route/CutBase.hpp"
#include <float.h>

bool SE_cut::validate_subset(const std::set<int> &S1, rut::Problem *prob) {
    if ( S1.size() <= 1 ||
        generatedSets.find(S1) != generatedSets.end() ) {
        return false;
    }
    std::set<int> setDiffRR;
    std::set_difference(prob->R.begin(), prob->R.end(), S1.begin(), S1.end(), std::inserter(setDiffRR, setDiffRR.end()));
    if (setDiffRR.size() == 0) {
        std::set<int> vistedD;
        std::set_intersection(prob->D.begin(), prob->D.end(), S1.begin(), S1.end(), std::inserter(vistedD, vistedD.end()));
        bool isAllPsVisited = true;
        for (int k: prob->K) {
            if (vistedD.find(prob->n_k[k]) == vistedD.end()) {
                continue;
            }
            if (S1.find(prob->h_k[k]) == S1.end()) {
                isAllPsVisited = false;
                break;
            }
        }
        return !isAllPsVisited;
    } else {
        return true;
    }
}

bool SE_cut::validate_LHS(const std::set<int> &S1, CutComposer *cc, const IloCplex::Callback::Context &context) {
    double lhs = 0.0;
    for (int i: S1) {
        for (int j: S1) {
            lhs += context.getRelaxationPoint(cc->x_ij[i][j]);
        }
    }
    lhs -= ((int) S1.size() - 1);
    if (lhs > 0) {
        return true;
    } else{
        return false;
    }
}

std::set<std::set<int>> SE_cut::solve_separationProb(CutComposer *cc, const IloCplex::Callback::Context &context) {
    rut::Problem *prob = cc->prob;
    std::vector<int> origins;
    for (int n0: prob->R) {
        if (n0 == prob->o || n0 == prob->d) {
            continue;
        }
        origins.push_back(n0);
    }
    std::set<std::set<int>> outputSets = get_maxCycles(origins, cc, context);
    //
    std::set<std::set<int>> validSets;
    for (std::set<int> S1: outputSets) {
        if (validate_subset(S1, prob) && validate_LHS(S1, cc, context)) {
            validSets.insert(S1);
        }
    }
    return validSets;
}

void SE_cut::set_LHS_Expr(IloExpr& lhs_expr, IloNumVar** x_ij, const std::set<int> &S1) {
    for (int i: S1) {
        for (int j: S1) {
            lhs_expr += x_ij[i][j];
        }
    }
    lhs_expr -= ((int) S1.size() - 1);
}

void SE_cut::add_cnsts2Model(const std::set<std::set<int>> &validSets,  CutComposer *cc, const IloCplex::Callback::Context &context) {
    for (std::set<int> S1: validSets) {
        generatedSets.insert(S1);
        IloExpr lhs_expr(cc->env);
        set_LHS_Expr(lhs_expr, cc->x_ij, S1);
        addUserCutwCust(context, lhs_expr);
        lhs_expr.end();
        cc->numGenCuts += 1;
    }
}

void SE_cut::add_cut(CutComposer *cc, const IloCplex::Callback::Context &context) {
    std::set<std::set<int>> validSets = solve_separationProb(cc, context);
    add_cnsts2Model(validSets, cc, context);
}

std::string SE_cut::add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context) {
    std::set<std::set<int>> validSets = solve_separationProb(cc, context);
    add_cnsts2Model(validSets, cc, context);
    //
    std::string addedCuts;
    for (std::set<int> S1: validSets) {
        addedCuts += "(";
        for (int i: S1) {
            addedCuts += std::to_string(i) + "-";
        }
        addedCuts += ");";
    }
    return addedCuts;
}

IloRangeArray SE_cut::get_detectedCuts(CutComposer *cc) {
    char buf[2048];
    IloRangeArray cnsts(cc->env);
    IloExpr lhs_expr(cc->env);
    for (std::set<int> S1: generatedSets) {
        lhs_expr.clear();
        sprintf(buf, "%s(%d)", cut_name.c_str(), cutCounter++);
        set_LHS_Expr(lhs_expr, cc->x_ij, S1);
        cnsts.add(lhs_expr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    generatedSets.clear();
    lhs_expr.end();
    return cnsts;
}



//
// The below codes are for dealing with cuts of the root node
//

std::set<std::set<int>> solve_maxSE_constructCycle_RC(double **x_ij, rut::Problem *prob) {
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
    int n1;
    const int n2 = prob->d;
    std::vector<int> path;
    for (int n0: prob->R) {
        if (n0 == prob->o || n0 == prob->d) {
            continue;
        }
        std::memset(visited, false, sizeof(visited));
        isCycle = false;
        //
        path.clear();
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
    }
    return outputSets;
}

IloRangeArray SE_cut::get_cut_cnsts(double **x_ij, CutComposer *cc) {
    rut::Problem *prob = cc->prob;
    std::set<std::set<int>> validSets;
    std::set<std::set<int>> outputSets = solve_maxSE_constructCycle_RC(x_ij, prob);
    for (std::set<int> S1: outputSets) {
        if (validate_subset(S1, prob)) {
            double lhs = 0.0;
            for (int i: S1) {
                for (int j: S1) {
                    lhs += x_ij[i][j];
                }
            }
            lhs -= ((int) S1.size() - 1);
            if (lhs > 0) {
                validSets.insert(S1);
            }
        }
    }
    //
    char buf[2048];
    IloRangeArray cnsts(cc->env);
    for (std::set<int> S1: validSets) {
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
