//
//  BaseMM.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 27/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef BaseMM_hpp
#define BaseMM_hpp

#include <cfloat>

#include <ilcplex/ilocplex.h>

#include "../Problem.hpp"
#include "../Etc.hpp"

#include "../Cuts/Base.hpp"


#define DEFAULT_BUFFER_SIZE 2048


class BaseMM {
public:
    Problem* prob;
    std::string logPath;
    IloEnv env;
    IloCplex* cplex;
    IloModel* cplexModel;
    IloNumVar** x_ij;
    IloNumVar* u_i;
    //
    BaseMM(Problem *prob, std::string logPath, char xType, bool isTightenModel);
    ~BaseMM();
    void get_x_ij(double** _x_ij);
    void get_u_i(double* _u_i);
    void start_fromGHSol(double** _x_ij, double* _u_i);
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
    void def_dvs(char xType);
    void def_FC_cnsts();
    void def_AT_cnsts();
    void def_CP_cnsts();
    void def_objF();
    void def_Ti_cnsts();
    void def_temp();
};


class ILP : public BaseMM {
public:
    ILP(Problem *prob, std::string logPath, bool isTightenModel);
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

class BnC : public BaseMM {
public:
    CutComposer *cc;
    //
    BnC(Problem *prob, std::string logPath, std::vector<std::string> cut_names, bool isTightenModel, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd, TimeTracker* tt);
    ~BnC() {
        delete cc;
    }
    int getNumGenCuts();
    double getTime4Sep();
    double getTime4FDV();
    int getNum4Sep();
};

class LP : public BaseMM {
public:
    std::ofstream* logStream;
    LP(Problem *prob, std::string logPath, bool isTightenModel);
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

class RC : public BaseMM {
public:
    CutComposer *cc;
    //
    RC(Problem *prob, std::string logPath, std::vector<std::string> cut_names);
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

#endif /* BaseMM_hpp */
