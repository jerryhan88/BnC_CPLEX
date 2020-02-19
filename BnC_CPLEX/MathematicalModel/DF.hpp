//
//  DF.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 20/1/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef DF_hpp
#define DF_hpp

#include <cfloat>

#include <ilcplex/ilocplex.h>

#include "../Problem.hpp"
#include "../Etc.hpp"

#define BUFFER_SIZE 2048   // same as DEFAULT_BUFFER_SIZE in BaseMM.hpp

class DF {
public:
    Problem* prob;
    std::string logPath;
    IloEnv env;
    IloCplex* cplex;
    IloModel* cplexModel;
    IloNumVar** x_ij;
    IloNumVar* u_i;
    //
    DF(Problem *prob, std::string logPath);
    ~DF();
    void get_x_ij(double** _x_ij);
    void get_x_ij(IloArray<IloNumArray>& _x_ij);
    void get_u_i(double* _u_i);
protected:
    void build_baseModel();
private:
    void def_dvs();
    void def_FC_cnsts();
    void def_AT_cnsts();
    void def_objF();
};

#endif /* DF_hpp */
