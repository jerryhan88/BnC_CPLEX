//
//  BaseMM.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 27/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "BaseMM.hpp"



void BaseMM::build_baseModel() {
    def_FC_cnsts();
    def_AT_cnsts();
    def_CP_cnsts();
    def_objF();
//    def_temp();
}

BaseMM::BaseMM(Problem *prob, std::string logPath, char xType, bool isTightenModel) {
    this->prob = prob;
    this->logPath = logPath;
    x_ij = new IloNumVar*[(*prob).N.size()];
    for (int i: (*prob).N) {
        x_ij[i] = new IloNumVar[(*prob).N.size()];
    }
    u_i = new IloNumVar[(*prob).N.size()];
    cplexModel = new IloModel(env);
    //
    def_dvs(xType);
    build_baseModel();
    if (isTightenModel) {
        def_Ti_cnsts();
    }
    //
    cplex = new IloCplex(*cplexModel);
}

BaseMM::~BaseMM() {
    for (int i = 0; i < (*prob).N.size(); i++)
        delete [] x_ij[i];
    delete [] x_ij;
    delete [] u_i;
    delete cplexModel;
    env.end();
}

void BaseMM::get_x_ij(double **_x_ij) {
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            _x_ij[i][j] = cplex->getValue(x_ij[i][j]);
        }
    }
}

void BaseMM::get_u_i(double *_u_i) {
    for (int i: (*prob).N) {
        _u_i[i] = cplex->getValue(u_i[i]);
    }
}


void BaseMM::start_fromGHSol(double** _x_ij, double* _u_i) {
    IloNumVarArray startVar(env);
    IloNumArray startVal(env);
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            startVar.add(x_ij[i][j]);
            startVal.add(_x_ij[i][j]);
        }
    }
    for (int i: (*prob).N) {
        startVar.add(u_i[i]);
        startVal.add(_u_i[i]);
    }
    IloCplex::MIPStartEffort effort = IloCplex::MIPStartSolveFixed;
    cplex->addMIPStart(startVar, startVal, effort);
    startVal.end();
    startVar.end();
}

std::vector<CutBase*> BaseMM::get_cutInstances(std::vector<std::string> cut_names, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd) {
    std::vector<CutBase*> cuts;
    for (std::string cut_name: cut_names) {
        std::string cutType = cut_name.substr(1);
        if (cutType == "SE") {
            cuts.push_back(new SE_cut(cut_name, cutManagerType, isLocalCutAdd));
        } else if (cutType == "CO") {
            assert(false);
//            cuts.push_back(new Cover_cut(cut_name, turnOnCutManager, isLocalCutAdd));
        } else if (cutType == "CA") {
            cuts.push_back(new CA_cut(cut_name, cutManagerType, isLocalCutAdd));
        } else if (cutType == "RS") {
            cuts.push_back(new RS_cut(cut_name, cutManagerType, isLocalCutAdd));
        } else if (cutType == "IP") {
            cuts.push_back(new IP_cut(cut_name, cutManagerType, isLocalCutAdd));
        } else {
            assert(false);
        }
    }
    return cuts;
}

void BaseMM::def_FC_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
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
    for (int i: (*prob).S) {
        if (i == (*prob).o || i == (*prob).d)
            continue;
        linExpr.clear();
        sprintf(buf, "iFS1(%d)", i);
        for (int j: (*prob).N) {
            if (i == j)
                continue;
            linExpr += x_ij[i][j];
        }
        cnsts.add(linExpr == 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
        //
        linExpr.clear();
        sprintf(buf, "iFS2(%d)", i);
        for (int j: (*prob).N) {
            if (i == j)
                continue;
            linExpr += x_ij[j][i];
        }
        cnsts.add(linExpr == 1);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    linExpr.clear();
    sprintf(buf, "DO");
    cnsts.add(x_ij[(*prob).d][(*prob).o] == 1);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.clear();
    sprintf(buf, "xF");
    for (int i: (*prob).N) {
        linExpr += x_ij[i][i];
    }
    cnsts.add(linExpr == 0);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    for (int k: (*prob).K) {
        linExpr.clear();
        sprintf(buf, "tFC(%d)", k);
        for (int j: (*prob).N) {
            linExpr += x_ij[(*prob).n_k[k]][j];
        }
        for (int j: (*prob).N) {
            linExpr -= x_ij[j][(*prob).h_k[k]];
        }
        cnsts.add(linExpr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    for (int i: (*prob).PD) {
        linExpr.clear();
        sprintf(buf, "FC_1(%d)", i);
        for (int j: (*prob).N) {
            linExpr += x_ij[i][j];
        }
        cnsts.add(linExpr <= 1);
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

void BaseMM::def_AT_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
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
        linExpr -= (*prob).M;
        for (int j: (*prob).N) {
            linExpr += (*prob).M * x_ij[(*prob).n_k[k]][j];
        }
        cnsts.add(linExpr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    for (int i: (*prob).S) {
        for (int j: (*prob).S) {
            linExpr.clear();
            sprintf(buf, "RR_P(%d,%d)", i, j);
            linExpr += (*prob).c_ij[i][j] * u_i[i];
            linExpr -= u_i[j];
            cnsts.add(linExpr <= 0);
            cnsts[cnsts.getSize() - 1].setName(buf);
        }
    }
    //
    linExpr.end();
    cplexModel->add(cnsts);
}

void BaseMM::def_CP_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    linExpr.clear();
    sprintf(buf, "bv");
    for (int k: (*prob).K) {
        for (int j: (*prob).N) {
            linExpr += (*prob).v_k[k] * x_ij[j][(*prob).n_k[k]];
        }
    }
    cnsts.add(linExpr <= (*prob).bv);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.clear();
    sprintf(buf, "bw");
    for (int k: (*prob).K) {
        for (int j: (*prob).N) {
            linExpr += (*prob).w_k[k] * x_ij[j][(*prob).n_k[k]];
        }
    }
    cnsts.add(linExpr <= (*prob).bw);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.clear();
    sprintf(buf, "DT");
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            if (i == (*prob).d && j == (*prob).o)
                continue;
            linExpr += (*prob).t_ij[i][j] * x_ij[i][j];
        }
    }
    cnsts.add(linExpr <= (*prob).bu);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.end();
    cplexModel->add(cnsts);
}

void BaseMM::def_objF() {
    IloExpr objF(env);
    for (int k: (*prob).K) {
        for (int j: (*prob).N) {
            objF += (*prob).r_k[k] * x_ij[j][(*prob).n_k[k]];
        }
    }
    cplexModel->add(IloMaximize(env, objF));
    objF.end();
}

void BaseMM::def_Ti_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    for (int i: (*prob).P) {
        linExpr.clear();
        sprintf(buf, "ED(%d)", i);
        for (int j: (*prob).N) {
            linExpr += x_ij[i][j];
        }
        for (int k: (*prob).K) {
            if ((*prob).h_k[k] == i) {
                for (int j: (*prob).N) {
                    linExpr -= x_ij[(*prob).n_k[k]][j];
                }
            }
        }
        cnsts.add(linExpr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    //
    linExpr.end();
    cplexModel->add(cnsts);
}

void BaseMM::def_temp() {
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    cnsts.add(x_ij[0][24] == 1);
    cnsts.add(x_ij[24][25] == 1);
    cnsts.add(x_ij[25][2] == 1);
    cnsts.add(x_ij[2][26] == 1);
    cnsts.add(x_ij[26][6] == 1);
    cnsts.add(x_ij[6][27] == 1);
    cnsts.add(x_ij[27][23] == 1);
    //
    linExpr.end();
    cplexModel->add(cnsts);
}

void BaseMM::def_dvs(char xType) {
    char buf[DEFAULT_BUFFER_SIZE];
    IloNumVar::Type ilo_var_type = xType == 'I' ? ILOINT : ILOFLOAT;
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            sprintf(buf, "x(%d)(%d)", i, j);
            x_ij[i][j] = IloNumVar(env, 0.0, 1.0, ilo_var_type, buf);
        }
    }
    for (int i: (*prob).N) {
        sprintf(buf, "u(%d)", i);
        u_i[i] = IloNumVar(env, 0.0, DBL_MAX, ILOFLOAT, buf);
    }
}

