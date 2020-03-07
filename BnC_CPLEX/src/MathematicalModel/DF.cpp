//
//  DF.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 20/1/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "DF.hpp"


void DF::build_baseModel() {
    def_FC_cnsts();
    def_AT_cnsts();
    def_objF();
}

DF::DF(rut::Problem *prob, std::string logPath) {
    this->prob = prob;
    this->logPath = logPath;
    x_ij = new IloNumVar*[(*prob).N.size()];
    for (int i: (*prob).N) {
        x_ij[i] = new IloNumVar[(*prob).N.size()];
    }
    u_i = new IloNumVar[(*prob).N.size()];
    cplexModel = new IloModel(env);
    //
    def_dvs();
    build_baseModel();
    //
    cplex = new IloCplex(*cplexModel);
}

DF::~DF() {
    for (int i = 0; i < (*prob).N.size(); i++)
        delete [] x_ij[i];
    delete [] x_ij;
    delete [] u_i;
    delete cplexModel;
    env.end();
}

void DF::get_x_ij(double **_x_ij) {
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            _x_ij[i][j] = cplex->getValue(x_ij[i][j]);
        }
    }
}

void DF::get_u_i(double *_u_i) {
    for (int i: (*prob).N) {
        _u_i[i] = cplex->getValue(u_i[i]);
    }
}

void DF::def_FC_cnsts() {
    char buf[BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    linExpr.clear();
    sprintf(buf, "iF1");
    for (int j: (*prob).N) {
        linExpr += x_ij[(*prob).o][j];
    }
    cnsts.add(linExpr == 1);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.clear();
    sprintf(buf, "iF2");
    for (int j: (*prob).N) {
        linExpr += x_ij[j][(*prob).d];
    }
    cnsts.add(linExpr == 1);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.clear();
    sprintf(buf, "DO");
    cnsts.add(x_ij[(*prob).d][(*prob).o] == 1);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    for (int i: (*prob).N) {
        sprintf(buf, "xF(%d)", i);
        cnsts.add(x_ij[i][i] == 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    for (int i: (*prob).PD) {
        linExpr.clear();
        sprintf(buf, "FC_1(%d)", i);
        for (int j: (*prob).N) {
            linExpr += x_ij[i][j];
        }
        cnsts.add(linExpr == 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
        //
        sprintf(buf, "FC(%d)", i);
        for (int j: (*prob).N) {
            linExpr -= x_ij[j][i];
        }
        cnsts.add(linExpr == 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    linExpr.end();
    cplexModel->add(cnsts);
}

void DF::def_AT_cnsts() {
    char buf[BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    sprintf(buf, "oA");
    cnsts.add(u_i[(*prob).o] == (*prob).al_i[(*prob).o]);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            if (i == (*prob).d && j == (*prob).o)
                continue;
            linExpr.clear();
            sprintf(buf, "AT(%d,%d)", i, j);
            linExpr += u_i[i] + (*prob).t_ij[i][j];
            linExpr -= u_i[j] + (*prob).M * (1 - x_ij[i][j]);
            cnsts.add(linExpr <= 0);
            cnsts[cnsts.getSize() - 1].setName(buf);
        }
    }
    //
    for (int i: (*prob).N) {
        sprintf(buf, "TW1(%d)", i);
        cnsts.add((*prob).al_i[i] <= u_i[i]);
        cnsts[cnsts.getSize() - 1].setName(buf);
        //
        sprintf(buf, "TW2(%d)", i);
        cnsts.add(u_i[i] <= (*prob).be_i[i]);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    for (int k: (*prob).K) {
        linExpr.clear();
        sprintf(buf, "WD_S(%d)", k);
        linExpr += u_i[(*prob).h_k[k]];
        linExpr -= u_i[(*prob).n_k[k]];
        cnsts.add(linExpr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    linExpr.clear();
    linExpr += u_i[(*prob).o];
    linExpr -= u_i[(*prob).d];
    cnsts.add(linExpr <= 0);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.end();
    cplexModel->add(cnsts);
}

void DF::def_objF() {
    IloExpr objF(env);
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            if (i == (*prob).d && j == (*prob).o)
                continue;
            objF += (*prob).t_ij[i][j] * x_ij[i][j];
        }
    }
    cplexModel->add(IloMinimize(env, objF));
    objF.end();
}

void DF::def_dvs() {
    char buf[BUFFER_SIZE];
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            sprintf(buf, "x(%d)(%d)", i, j);
            x_ij[i][j] = IloNumVar(env, 0.0, 1.0, ILOINT, buf);
        }
    }
    for (int i: (*prob).N) {
        sprintf(buf, "u(%d)", i);
        u_i[i] = IloNumVar(env, 0.0, DBL_MAX, ILOFLOAT, buf);
    }
}
