//
//  Capacity.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 16/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Base.hpp"
#include "../NetworkFlow/FordFulkersonAlgo.hpp"

int get_numP(std::set<int> S1, Problem* prob) {
    int numP = 0;
    for (int i: S1) {
        if (prob->P.find(i) != prob->P.end()) {
            numP += 1;
        }
    }
    return numP;
}

int get_numVisitVol(std::set<int> S1, Problem* prob) {
    int numP = get_numP(S1, prob);
    int numVisitVol = std::floor( numP / (float) prob->bv);
    numVisitVol = numVisitVol > 1 ? numVisitVol : 1;
    return numVisitVol;
}

int get_numVisitWei(std::set<int> S1, Problem* prob) {
    int numP = get_numP(S1, prob);
    int numVisitWei = std::floor( numP / (float) prob->bw);
    numVisitWei = numVisitWei > 1 ? numVisitWei : 1;
    return numVisitWei;
}

bool Capacity_cut::valid_subset(const std::set<int>& S1, Problem *prob) {
    if ( S1.size() <= 1 ||
        generatedSets.find(S1) != generatedSets.end() ) {
        return false;
    }
    if ( S1.find(prob->o) != S1.end() ||
        S1.find(prob->d) != S1.end() ) {
        return false;
    }
    return true;
}

std::set<std::set<int>> Capacity_cut::solve_separationProb(double** x_ij, CutComposer* cc) {
    Problem *prob = cc->prob;
    std::set<std::set<int>> validSets;
    std::set<std::set<int>> outputSets = solve_multipleMaxFlow(x_ij, prob);
    for (std::set<int> S1: outputSets) {
        if (valid_subset(S1, prob)) {
            double lhs = 0.0;
            for (int i: S1) {
                for (int j: S1) {
                    lhs += x_ij[i][j];
                }
            }
            int numVisitVol, numVisitWei, maxCapa;
            numVisitVol = get_numVisitVol(S1, prob);
            numVisitWei = get_numVisitWei(S1, prob);
            maxCapa = numVisitVol > numVisitWei ? numVisitVol : numVisitWei;
            lhs -= ((int) S1.size() - maxCapa);
            if (lhs > 0) {
                validSets.insert(S1);
            }
        }
    }
    return validSets;
}

void Capacity_cut::set_LHS_Expr(IloExpr& lhs_expr, IloNumVar** x_ij, const std::set<int>& S1, Problem* prob) {
    for (int i: S1) {
        for (int j: S1) {
            lhs_expr += x_ij[i][j];
        }
    }
    int numVisitVol, numVisitWei, maxCapa;
    numVisitVol = get_numVisitVol(S1, prob);
    numVisitWei = get_numVisitWei(S1, prob);
    maxCapa = numVisitVol > numVisitWei ? numVisitVol : numVisitWei;
    lhs_expr -= ((int) S1.size() - maxCapa);
}

IloRangeArray Capacity_cut::get_cut_cnsts(double** x_ij, CutComposer* cc) {
    std::set<std::set<int>> validSets = solve_separationProb(x_ij, cc);
    char buf[2048];
    IloRangeArray cnsts(cc->env);
    for (std::set<int> S1: validSets) {
        sprintf(buf, "%s(%d)", cut_name.c_str(), (int) generatedSets.size());
        generatedSets.insert(S1);
        IloExpr lhs_expr(cc->env);
        set_LHS_Expr(lhs_expr, cc->x_ij, S1, cc->prob);
        cnsts.add(lhs_expr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
        lhs_expr.end();
    }
    return cnsts;
}

void Capacity_cut::add_cnsts2Model(const std::set<std::set<int>>& validSets,  CutComposer* cc, const IloCplex::Callback::Context& context) {
    for (std::set<int> S1: validSets) {
        generatedSets.insert(S1);
        IloExpr lhs_expr(cc->env);
        set_LHS_Expr(lhs_expr, cc->x_ij, S1, cc->prob);
        context.addUserCut(lhs_expr <= 0, IloCplex::UseCutForce, IloFalse);
        lhs_expr.end();
    }
}

void Capacity_cut::add_cut(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context) {
    std::set<std::set<int>> validSets = solve_separationProb(x_ij, cc);
    add_cnsts2Model(validSets, cc, context);
}

std::string Capacity_cut::add_cut_wLogging(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context) {
    std::set<std::set<int>> validSets = solve_separationProb(x_ij, cc);
    add_cnsts2Model(validSets, cc, context);
    std::string addedCuts;
    for (std::set<int> S1: validSets) {
        int numVisitVol, numVisitWei, maxCapa;
        numVisitVol = get_numVisitVol(S1, cc->prob);
        numVisitWei = get_numVisitWei(S1, cc->prob);
        maxCapa = numVisitVol > numVisitWei ? numVisitVol : numVisitWei;
        addedCuts += "(";
        for (int i: S1) {
            addedCuts += std::to_string(i) + "-";
        }
        addedCuts += ")" + std::to_string(maxCapa)+ ";";
    }
    //
    return addedCuts;
}




