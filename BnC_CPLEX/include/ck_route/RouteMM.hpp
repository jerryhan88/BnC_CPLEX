//
//  RouteMM.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 1/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef RouteMM_h
#define RouteMM_h


#include <cfloat>

#include <ilcplex/ilocplex.h>

#include "Problem.hpp"
#include "CutBase.hpp"

#include "ck_util/util.hpp"

#define DEFAULT_BUFFER_SIZE 2048

namespace rmm {

class RouteMM {
public:
    rut::Problem *prob;
    std::string logPath;
    IloEnv env;
    IloCplex *cplex;
    IloModel *cplexModel;
    IloNumVar **x_ij, *u_i;
    //
    RouteMM(rut::Problem *prob, std::string logPath, char vType, bool isTightenModel) {
        this->prob = prob;
        this->logPath = logPath;
        //
        x_ij = RouteMM::gen_x_ij(prob, env, vType);
        u_i = RouteMM::gen_u_i(prob, env);
        
        cplexModel = new IloModel(env);
        build_baseModel();
        if (isTightenModel) {
            def_Ti_cnsts();
        }
        //
        cplex = new IloCplex(*cplexModel);
    }
    ~RouteMM() {
        for (int i = 0; i < prob->N.size(); i++)
            delete [] x_ij[i];
        delete [] x_ij; delete [] u_i;
        delete cplexModel; delete cplex;
        env.end();
    }
    static IloNumVar** gen_x_ij(rut::Problem *prob, IloEnv &env, char vType) {
        IloNumVar::Type ilo_vType = vType == 'I' ? ILOINT : ILOFLOAT;
        //
        char buf[DEFAULT_BUFFER_SIZE];
        IloNumVar **x_ij = new IloNumVar*[prob->N.size()];
        for (int i: prob->N) {
            x_ij[i] = new IloNumVar[prob->N.size()];
            for (int j: prob->N) {
                sprintf(buf, "x(%d)(%d)", i, j);
                x_ij[i][j] = IloNumVar(env, 0.0, 1.0, ilo_vType, buf);
            }
        }
        return x_ij;
    }
    static IloNumVar* gen_u_i(rut::Problem *prob, IloEnv &env) {
        char buf[DEFAULT_BUFFER_SIZE];
        IloNumVar *u_i = new IloNumVar[prob->N.size()];
        for (int i: prob->N) {
            sprintf(buf, "u(%d)", i);
            u_i[i] = IloNumVar(env, 0.0, DBL_MAX, ILOFLOAT, buf);
        }
        return u_i;
    }
    void get_x_ij(double** _x_ij);
    void get_x_ij(IloArray<IloNumArray>& _x_ij);
    void get_u_i(double* _u_i);
    void set_initSol(double** _x_ij, double* _u_i);
    virtual int getNumGenCuts() {
           throw "Should override getNumGenCuts()";
    }
    virtual double getTime4Sep() {
        throw "Should override getTime4Sep()";
    }
    virtual double getTime4FDV() {
        throw "Should override getTime4FDV()";
    }
    virtual int getNum4Sep() {
        throw "Should override getNum4Sep()";
    }
protected:
    void build_baseModel();
    std::vector<CutBase*> get_cutInstances(std::vector<std::string> cut_names, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd);
private:
    void def_FC_cnsts();
    void def_AT_cnsts();
    void def_CP_cnsts();
    void def_objF();
    void def_Ti_cnsts();
    void def_temp();
};


class ILP : public RouteMM {
public:
    ILP(rut::Problem *prob, std::string logPath, bool isTightenModel);
    int getNumGenCuts() {
        return -1;
    }
    double getTime4Sep() {
        return -1.0;
    }
    double getTime4FDV() {
        return -1.0;
    }
    int getNum4Sep() {
        return -1;
    }
};

class BnC : public RouteMM {
public:
    CutComposer *cc;
    //
    BnC(rut::Problem *prob, std::string logPath, std::vector<std::string> cut_names, bool isTightenModel, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd, TimeTracker* tt);
    ~BnC() {
        delete cc;
    }
    int getNumGenCuts();
    double getTime4Sep();
    double getTime4FDV();
    int getNum4Sep();
    void add_detectedCuts2MM();
};

class LP : public RouteMM {
public:
    std::ofstream* logStream;
    LP(rut::Problem *prob, std::string logPath, bool isTightenModel);
    ~LP();
    int getNumGenCuts() {
        return -1;
    }
    double getTime4Sep() {
        return -1.0;
    }
    double getTime4FDV() {
        return -1.0;
    }
    int getNum4Sep() {
        return -1;
    }
};

class RC : public RouteMM {
public:
    CutComposer *cc;
    //
    RC(rut::Problem *prob, std::string logPath, std::vector<std::string> cut_names);
    //
    ~RC() {
        delete cc;
    }
    void solve();
    int getNumGenCuts() {
        return -1;
    }
    double getTime4FDV() {
        return -1.0;
    }
    double getTime4Sep() {
        return -1.0;
    }
    int getNum4Sep() {
        return -1;
    }
};

}


#endif /* RouteMM_h */
