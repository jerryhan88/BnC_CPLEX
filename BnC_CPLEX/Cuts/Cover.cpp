//
//  Cover.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 14/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Base.hpp"

double min(double a, double b) { return (a < b)? a : b; }

void knapSack(const std::vector<double>& values, const std::vector<int>& weights, int minWeight, double* obj, bool excludedEdge[]) {
    int sumWeights = 0, numItems = (int) values.size();
    double maxValue = 0.0;
    for (int i = 0; i < numItems; i++) {
        maxValue += values[i];
        sumWeights += weights[i];
    }
    double** K = new double*[numItems + 1];
    for (int i = 0; i < numItems + 1; i++) {
        K[i] = new double[sumWeights + 1];
    }
    for (int i = 0; i <= numItems; i++) {
        for (int w = 0; w <= sumWeights; w++) {
            if (i == 0 || w == sumWeights) {
                K[i][w] = maxValue;
            } else if ( (sumWeights - weights[i - 1]) >= w ) {
                K[i][w] = min(K[i - 1][w + weights[i - 1]] - values[i - 1],  K[i - 1][w]);
            } else {
                K[i][w] = K[i - 1][w];
            }
        }
    }
    *obj = K[numItems][minWeight];
    //
    int n = numItems, M = minWeight;
    while (n !=0) {
        if (K[n][M] != K[n - 1][M]) {
            excludedEdge[n - 1] = true;
            M += weights[n - 1];
        }
        n--;
    }
    for (int i = 0; i < numItems + 1; i++) {
        delete [] K[i];
    }
    delete [] K;
}

bool Cover_cut::valid_subset(const std::set<edge>& S1) {
    if (generatedSets.find(S1) != generatedSets.end() ) {
        return false;
    } else {
        return true;
    }
}

std::set<edge> Cover_cut::solve_separationProb(double** x_ij, CutComposer* cc) {
    Problem *prob = cc->prob;
    std::set<edge> selectedEdges;
    double obj;
    int coveredWeight = 0;
    std::vector<edge> targetedEdges;
    std::vector<double> values;
    std::vector<int> weights;
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            if (x_ij[i][j] == 1.0) {
                selectedEdges.insert(edge {i, j});
                coveredWeight += (*prob).t_ij[i][j];
            } else {
                if (x_ij[i][j] == 0.0) {
                    continue;
                }
                targetedEdges.push_back(edge {i, j});
                values.push_back(1 - x_ij[i][j]);
                weights.push_back((*prob).t_ij[i][j]);
            }
        }
    }
    int minWeight = (*prob).bu - coveredWeight + 1;
    bool excludedEdge[values.size()];
    std::memset(excludedEdge, false, sizeof(excludedEdge));
    knapSack(values, weights, minWeight, &obj, excludedEdge);
    for (int i = 0; i < values.size(); i++) {
        if (!excludedEdge[i]) {
            selectedEdges.insert(targetedEdges[i]);
        }
    }
    if (obj < 1.0) {
        if (!valid_subset(selectedEdges)) {
            selectedEdges.clear();
        }
    } else {
        selectedEdges.clear();
    }
    return selectedEdges;
}

void Cover_cut::set_LHS_Expr(IloExpr& lhs_expr, IloNumVar** x_ij, const std::set<edge>& S1) {
    for (edge e: S1) {
        lhs_expr += x_ij[e.first][e.second];
    }
    lhs_expr -= ((int) S1.size() - 1);
}

IloRangeArray Cover_cut::get_cut_cnsts(double** x_ij, CutComposer* cc) {
    std::set<edge> S1 = solve_separationProb(x_ij, cc);
    IloRangeArray cnsts(cc->env);
    if (S1.size() != 0) {
        char buf[2048];
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

void Cover_cut::add_cnsts2Model(const std::set<edge>& S1, CutComposer* cc, const IloCplex::Callback::Context& context) {
    generatedSets.insert(S1);
    IloExpr lhs_expr(cc->env);
    set_LHS_Expr(lhs_expr, cc->x_ij, S1);
    context.addUserCut(lhs_expr <= 0, IloCplex::UseCutForce, IloFalse);
    lhs_expr.end();
}

void Cover_cut::add_cut(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context) {
    std::set<edge> S1 = solve_separationProb(x_ij, cc);
    if (S1.size() != 0) {
        add_cnsts2Model(S1, cc, context);
    }
}

std::string Cover_cut::add_cut_wLogging(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context) {
    std::set<edge> S1 = solve_separationProb(x_ij, cc);
    std::string addedCuts;
    if (S1.size() != 0) {
        add_cnsts2Model(S1, cc, context);
        addedCuts += "(";
        for (edge e: S1) {
            addedCuts += std::to_string(e.first) + "#" + std::to_string(e.second) + "-";
        }
        addedCuts += ");";
    }
    return addedCuts;
}
