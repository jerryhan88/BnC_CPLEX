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

#include "../../include/ck_route/Other.hpp"

#include "ck_util/util.hpp"

#define BUFFER_SIZE 2048   // same as DEFAULT_BUFFER_SIZE in BaseMM.hpp


class DF {
public:
    rut::Problem* prob;
    TimeTracker* tt;
    unsigned long time_limit_sec;
    std::string logPath;
    std::string lpPath;
    IloEnv env;
    IloCplex* cplex;
    IloModel* cplexModel;
    IloNumVar** x_ij;
    IloNumVar* u_i;
    //
    DF(rut::Problem *prob, TimeTracker* tt, unsigned long time_limit_sec, std::string logPath, std::string lpPath) {
        this->prob = prob;
        this->tt = tt;
        this->time_limit_sec = time_limit_sec;
        this->logPath = logPath;
        this->lpPath = lpPath;
        //
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
    ~DF() {
        for (int i = 0; i < (*prob).N.size(); i++)
            delete [] x_ij[i];
        delete [] x_ij;
        delete [] u_i;
        delete cplexModel;
        env.end();
    }
    //
    rut::Solution* solve();
    //
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
