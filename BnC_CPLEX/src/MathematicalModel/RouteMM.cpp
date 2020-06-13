//
//  RouteMM.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 1/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ck_route/Router.hpp"

rut::Solution* rmm::RouteMM::solve() {
    run();
    //
    rut::Solution *sol = new rut::Solution(prob);
    //
    try {
        sol->objV = cplex->getObjValue();
        sol->gap = cplex->getMIPRelativeGap();
        sol->cpuT = tt->get_elapsedTimeCPU();
        sol->wallT = tt->get_elapsedTimeWall();
        //
        char note[2048];
        sprintf(note,"\"{\'numNodes\': %lld, \'numGenCuts\': %d, \'time4Sep\': %f,\'time4FDV\': %f, \'num4Sep\': %d}\"",
                        cplex->getNnodes64(),
                        getNumGenCuts(),
                        getTime4Sep(),
                        getTime4FDV(),
                        getNum4Sep());
        sol->note = std::string(note);
        //
        for (int i: prob->N) {
            for (int j: prob->N) {
                sol->x_ij[i][j] = cplex->getValue(x_ij[i][j]);
            }
            sol->u_i[i] = cplex->getValue(u_i[i]);
        }
    } catch (IloCplex::Exception e) {
        std::cout << "no incumbent until the time limit" << std::endl;
        std::fstream fout;
        fout.open(lpPath, std::ios::out);
        fout << "\t No incumbent until the time limit, " << time_limit_sec << "seconds" << "\n";
        fout.close();
    }
    return sol;
}

void rmm::RouteMM::run() {
    try {
        cplex->solve();
    } catch(IloException& e) {
       std::cerr << "Concert exception caught" << e << std::endl;
       throw;
    }
    if (cplex->getStatus() == IloAlgorithm::Infeasible) {
        cplex->exportModel(lpPath.c_str());
        throw "Infeasible";
    }
}

void rmm::RouteMM::get_x_ij(double **_x_ij) {
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            _x_ij[i][j] = cplex->getValue(x_ij[i][j]);
        }
    }
}

void rmm::RouteMM::get_u_i(double *_u_i) {
    for (int i: (*prob).N) {
        _u_i[i] = cplex->getValue(u_i[i]);
    }
}

void rmm::RouteMM::build_baseModel() {
    def_preprocessing();
    //
    def_FC_cnsts();
    def_AT_cnsts();
    def_CP_cnsts();
    def_objF();
//    def_temp();
}

void rmm::RouteMM::def_preprocessing() {
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    //
    std::set<iiTup> invalid_dvX;
    for (int i0 = 0; i0 < prob->_R.size() - 1; i0++) {
        int n0 = prob->_R[i0];
        double al_n0 = prob->al_i[n0];
        int n2 = prob->_R[i0 + 1];
        double be_n2 = prob->be_i[n2];
        for (int n1: prob->P) {
            if (al_n0 + prob->t_ij[n0][n1] + prob->t_ij[n1][n2] > be_n2) {
                invalid_dvX.insert(std::make_tuple(n0, n1));
                invalid_dvX.insert(std::make_tuple(n1, n2));
            }
        }
        for (int n1: prob->D) {
            if (al_n0 + prob->t_ij[n0][n1] + prob->t_ij[n1][n2] > be_n2) {
                invalid_dvX.insert(std::make_tuple(n0, n1));
                invalid_dvX.insert(std::make_tuple(n1, n2));
            }
        }
    }
    for (int i: prob->N) {
        for (int j: prob->N) {
            if (i == prob->d && j == prob->o) {
                continue;
            }
            if (prob->al_i[i] + prob->t_ij[i][j] > prob->be_i[j]) {
                invalid_dvX.insert(std::make_tuple(i, j));
            }
        }
    }
    int o = prob->o;
    for (int i: prob->N) {
        if (i == prob->d) {
            continue;
        }
        invalid_dvX.insert(std::make_tuple(i, o));
    }
    int d = prob->d;
    for (int j: prob->N) {
        if (j == prob->o) {
            continue;
        }
        invalid_dvX.insert(std::make_tuple(d, j));
    }
    for (int k : prob->K) {
        int i = prob->n_k[k];
        int j = prob->h_k[k];
        invalid_dvX.insert(std::make_tuple(i, j));
    }
    int i, j;
    for (iiTup ele: invalid_dvX) {
        std::tie(i, j) = ele;
        cnsts.add(x_ij[i][i] == 0);
        sprintf(buf, "pP(%d)(%d)", i, j);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
    cplexModel->add(cnsts);
}

void rmm::RouteMM::set_initSol(double** _x_ij, double* _u_i) {
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

std::vector<CutBase*> rmm::RouteMM::get_cutInstances(std::vector<std::string> cut_names, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd) {
    std::vector<CutBase*> cuts;
    if (cut_names.size() == 4) {
        cuts.push_back(new ALL_cut("ALL", cutManagerType, isLocalCutAdd));
    } else {
        for (std::string cut_name: cut_names) {
            std::string cutType = cut_name.substr(1);
            if (cutType == "SE") {
                cuts.push_back(new SE_cut(cut_name, cutManagerType, isLocalCutAdd));
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
    }
    return cuts;
}

void rmm::RouteMM::def_FC_cnsts() {
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
    for (int i: (*prob).R) {
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
    sprintf(buf, "CF");  // Circular Flow for the branch-and-cut algorithm
    cnsts.add(x_ij[(*prob).d][(*prob).o] == 1);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    for (int i: (*prob).N) {
        sprintf(buf, "NS(%d)", i);  // No Self Flow; tightening bounds
        cnsts.add(x_ij[i][i] == 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
    }
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

void rmm::RouteMM::def_AT_cnsts() {
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
    for (int i: (*prob).R) {
        for (int j: (*prob).R) {
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

void rmm::RouteMM::def_CP_cnsts() {
    char buf[DEFAULT_BUFFER_SIZE];
    IloRangeArray cnsts(env);
    IloExpr linExpr(env);
    //
    linExpr.clear();
    sprintf(buf, "VL");  // Volume Limit
    for (int k: (*prob).K) {
        for (int j: (*prob).N) {
            linExpr += (*prob).v_k[k] * x_ij[j][(*prob).n_k[k]];
        }
    }
    cnsts.add(linExpr <= (*prob).bv);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.clear();
    sprintf(buf, "WL");  // Weight Limit
    for (int k: (*prob).K) {
        for (int j: (*prob).N) {
            linExpr += (*prob).w_k[k] * x_ij[j][(*prob).n_k[k]];
        }
    }
    cnsts.add(linExpr <= (*prob).bw);
    cnsts[cnsts.getSize() - 1].setName(buf);
    //
    linExpr.clear();
    sprintf(buf, "TL");  // Time Limit
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

void rmm::RouteMM::def_objF() {
    IloExpr objF(env);
    for (int k: (*prob).K) {
        for (int j: (*prob).N) {
            objF += (*prob).r_k[k] * x_ij[j][(*prob).n_k[k]];
        }
    }
    cplexModel->add(IloMaximize(env, objF));
    objF.end();
}

void rmm::RouteMM::def_Ti_cnsts() {
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

void rmm::RouteMM::def_temp() {
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
