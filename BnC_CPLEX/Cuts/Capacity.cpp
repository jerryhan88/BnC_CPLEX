//
//  Capacity.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 16/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Base.hpp"


int get_numP(std::set<int> S1, Problem* prob) {
    int numP = 0;
    for (int i: S1) {
        if (prob->P.find(i) != prob->P.end()) {
            numP += 1;
        }
    }
    return numP;
}

double get_vS(std::set<int> S1, Problem* prob) {
    double vS = 0.0;
    for (int n0: S1) {
        vS += (*prob).v_i[n0];
    }
    return vS;
}

double get_wS(std::set<int> S1, Problem* prob) {
    double wS = 0.0;
    for (int n0: S1) {
        wS += (*prob).w_i[n0];
    }
    return wS;
}

double get_LHS_CA(const std::set<int>& S, double** x_ij, Problem* prob) {
    double lhs = 0.0;
    for (int i: S) {
        for (int j: S) {
            lhs += x_ij[i][j];
        }
    }
    int _vS, _wS, maxCapa;
    _vS = std::floor(get_vS(S, prob) / (*prob).bv);
    _wS = std::floor(get_wS(S, prob) / (*prob).bw);
    maxCapa = _vS > _wS ? _vS : _wS;
    lhs -= ((int) S.size() - maxCapa);
    return lhs;
}

std::set<std::set<int>> solve_maxCA(double** x_ij, Problem* prob) {
    std::set<std::set<int>> outputSets;
    run_TB4CG(x_ij, prob, get_LHS_CA, outputSets);
    return outputSets;
}

Capacity_cut::Capacity_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd) : CutBase(ch_name, cutManagerType, isLocalCutAdd) {
    if (ch_name.substr(0, 1) == "e") {
        separationAlgo = solve_multipleMaxFlow;
    } else {
        assert(ch_name.substr(0, 1) == "h");
        separationAlgo = solve_maxCA;
    }
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
    std::set<std::set<int>> outputSets = separationAlgo(x_ij, prob);
    for (std::set<int> S1: outputSets) {
        if (valid_subset(S1, prob)) {
            double lhs = get_LHS_CA(S1, x_ij, prob);
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
    int _vS, _wS, maxCapa;
    _vS = std::floor(get_vS(S1, prob) / (*prob).bv);
    _wS = std::floor(get_wS(S1, prob) / (*prob).bw);
    maxCapa = _vS > _wS ? _vS : _wS;
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
        addUserCutwCust(context, lhs_expr);
        lhs_expr.end();
        cc->numGenCuts += 1;
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
        numVisitVol = get_vS(S1, cc->prob);
        numVisitWei = get_wS(S1, cc->prob);
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




