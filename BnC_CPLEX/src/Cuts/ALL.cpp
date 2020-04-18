//
//  ALL.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 10/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ck_route/CutBase.hpp"


void ALL_cut::add_cut(CutComposer *cc, const IloCplex::Callback::Context &context) {
    std::set<std::set<int>> validSets_se, validSets_ca, validSets_rs;
    //
    rut::Problem *prob = cc->prob;
    int numNodes = (int) prob->N.size();
    bool visited[numNodes], isCycle;
    const int n2 = prob->d;
    std::vector<int> path;
    for (int n0: prob->R) {
        if (n0 == prob->d) {
            continue;
        }
        isCycle = get_cyclicPath(n0, n2, visited, numNodes, path, cc, context);
        if (!isCycle) {
            continue;
        }
        std::set<int> S1(path.begin(), path.end());
        if (n0 == prob->o) {
            if (rs_cut->validate_subset(S1, prob) && rs_cut->validate_LHS(S1, prob, cc, context)) {
                validSets_rs.insert(S1);
            }
        } else {
            if (ca_cut->validate_subset(S1, prob) && ca_cut->validate_LHS(S1, prob, cc, context)) {
                validSets_ca.insert(S1);
            } else if (se_cut->validate_subset(S1, prob) && se_cut->validate_LHS(S1, cc, context)) {
                validSets_se.insert(S1);
            }
        }
    }
    se_cut->add_cnsts2Model(validSets_se, cc, context);
    ca_cut->add_cnsts2Model(validSets_ca, cc, context);
    rs_cut->add_cnsts2Model(validSets_rs, cc, context);
    //
    ip_cut->add_cut(cc, context);
}

IloRangeArray ALL_cut::get_detectedCuts(CutComposer *cc) {
    char buf[2048];
    IloRangeArray cnsts(cc->env);
    IloExpr lhs_expr(cc->env);
    //
    for (std::set<int> S1: se_cut->generatedSets) {
        lhs_expr.clear();
        sprintf(buf, "%s(%d)", se_cut->cut_name.c_str(), cutCounter++);
        se_cut->set_LHS_Expr(lhs_expr, cc->x_ij, S1);
        cnsts.add(lhs_expr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    se_cut->generatedSets.clear();
    //
    for (std::set<int> S1: ca_cut->generatedSets) {
        lhs_expr.clear();
        sprintf(buf, "%s(%d)", ca_cut->cut_name.c_str(), cutCounter++);
        ca_cut->set_LHS_Expr(lhs_expr, cc->x_ij, S1, cc->prob);
        cnsts.add(lhs_expr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    ca_cut->generatedSets.clear();
    //
    for (std::set<int> S1: rs_cut->generatedSets) {
        lhs_expr.clear();
        sprintf(buf, "%s(%d)", rs_cut->cut_name.c_str(), cutCounter++);
        rs_cut->set_LHS_Expr(lhs_expr, cc->x_ij, S1, cc->prob);
        cnsts.add(lhs_expr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    rs_cut->generatedSets.clear();
    //
    for (std::set<edge> S1: ip_cut->generatedSets) {
        lhs_expr.clear();
        sprintf(buf, "%s(%d)", ip_cut->cut_name.c_str(), cutCounter++);
        ip_cut->set_LHS_Expr(lhs_expr, cc->x_ij, S1);
        cnsts.add(lhs_expr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    ip_cut->generatedSets.clear();
    //
    lhs_expr.end();
    return cnsts;
}
