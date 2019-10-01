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
    void get_x_ij(double **_x_ij);
    void get_u_i(double *_u_i);
protected:
    void build_baseModel();
    std::vector<CutBase*> get_cutInstances(std::vector<std::string> cut_names);
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
};

class BnC : public BaseMM {
public:
    CutComposer *cc;
    //
    BnC(Problem *prob, std::string logPath, std::vector<std::string> cut_names);
    ~BnC() {
        delete cc;
    }
};

class LP : public BaseMM {
public:
    std::ofstream* logStream;
    LP(Problem *prob, std::string logPath, bool isTightenModel);
    ~LP();
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
};

#endif /* BaseMM_hpp */
