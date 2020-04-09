//
//  RouteMM.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 1/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef RouteMM_h
#define RouteMM_h

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <cfloat>

#include <ilcplex/ilocplex.h>

#include "Other.hpp"
#include "CutBase.hpp"
#include "Sequencer.h"

#include "ck_util/util.hpp"

#define DEFAULT_BUFFER_SIZE 2048

typedef std::pair< int, std::vector<int> > tid_bestSeq;

class Router {
public:
    rut::Problem* prob;
    TimeTracker* tt;
    unsigned long time_limit_sec;
    std::string logPath;
    IloCplex *cplex;
    //
    Router(rut::Problem *prob, TimeTracker* tt, unsigned long time_limit_sec, std::string logPath) {
        this->prob = prob;
        this->tt = tt;
        this->time_limit_sec = time_limit_sec;
        this->logPath = logPath;
    }
    virtual void run() {
           throw "Should override run()";
    }
    virtual rut::Solution* solve() {
           throw "Should override solve()";
    }
};

namespace rgh {

class InsertionHeuristic: public Router {
public:
    double **x_ij, *u_i;
    double limDet, limVol, limWei;
    double curReward, curVol, curWei;
    std::set<int> candiTasks, visitedWH;
    std::vector<int> partialSeq;
    int seqBeginIndex4Search;
    //
    InsertionHeuristic(rut::Problem *prob, TimeTracker* tt, unsigned long time_limit_sec, std::string logPath): Router(prob, tt, time_limit_sec, logPath)  {
        x_ij = new double*[prob->N.size()];
        u_i = new double[prob->N.size()];
        for (int i: prob->N) {
            x_ij[i] = new double[(*prob).N.size()];
        }
        limDet = (*prob).bu;
        limVol = (*prob).bv;
        limWei = (*prob).bw;
        initIH();
    }
    ~InsertionHeuristic() {
        for (int i = 0; i < prob->N.size(); i++)
            delete [] x_ij[i];
        delete [] x_ij;
        delete [] u_i;
    }
    void run();
    rut::Solution* solve();
    void initIH();
    tid_bestSeq findBestTaskSeq();
    double get_objV();
};
}

namespace rmm {

class RouteMM: public Router {
public:
    std::string lpPath;    //
    IloEnv env;
    IloModel *cplexModel;
    IloNumVar **x_ij, *u_i;
    //
    RouteMM(rut::Problem *prob, TimeTracker* tt, unsigned long time_limit_sec, std::string logPath, std::string lpPath, char vType, bool isTightenModel): Router(prob, tt, time_limit_sec, logPath) {
        this->lpPath = lpPath;
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
    void run();
    rut::Solution* solve();
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
    ILP(rut::Problem *prob, TimeTracker* tt, unsigned long time_limit_sec,
        std::string logPath, std::string lpPath, bool isTightenModel);
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
    BnC(rut::Problem *prob, TimeTracker* tt, unsigned long time_limit_sec,
        std::string logPath, std::string lpPath,
        std::vector<std::string> cut_names, bool isTightenModel,
        IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd) : RouteMM(prob, tt, time_limit_sec, logPath, lpPath, 'I', isTightenModel) {
        std::vector<CutBase*> cuts = get_cutInstances(cut_names, cutManagerType, isLocalCutAdd);
        cc = new CutComposer(prob, cuts, env, x_ij, logPath, tt);
        CPXLONG contextMask = 0;
        contextMask |= IloCplex::Callback::Context::Id::Relaxation;
        cplex->use(cc, contextMask);
        //
        cplex->setOut(env.getNullStream());
        if (logPath != "") {
            std::string _header("nCnt,bestObj,bestBound");
            for (CutBase *c: cc->cuts) {
                _header += "," + c->cut_name;
            }
            _header += ",note";
            char header[_header.size() + 1];
            std::strcpy(header, _header.c_str());
            createCSV(logPath, header);
        }
        
    }
    ~BnC() {
        delete cc;
    }
    int getNumGenCuts() {
        return cc->numGenCuts;
    }
    double getTime4Sep() {
        return cc->time4Sep;
    }
    double getTime4FDV() {
        return cc->time4FDV;
    }
    int getNum4Sep() {
        return cc->num4Sep;
    }
    void add_detectedCuts2MM() {
        for (CutBase *c: cc->cuts) {
            cplexModel->add(c->get_detectedCuts(cc));
        }
    }
};

class LP : public RouteMM {
public:
    std::ofstream* logStream;
    LP(rut::Problem *prob, TimeTracker* tt, unsigned long time_limit_sec, std::string logPath, std::string lpPath, bool isTightenModel): RouteMM(prob, tt, time_limit_sec, logPath, lpPath, 'L', isTightenModel) {
        logStream = new std::ofstream (logPath.c_str(), std::ios_base::app);
        cplex->setOut(*logStream);
    }
    ~LP() {
        delete logStream;
    }
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
    RC(rut::Problem *prob, TimeTracker* tt, unsigned long time_limit_sec, std::string logPath, std::string lpPath, std::vector<std::string> cut_names): RouteMM(prob, tt, time_limit_sec, logPath, lpPath, 'L', true)  {
        std::vector<CutBase*> cuts = get_cutInstances(cut_names, IloCplex::UseCutForce, IloFalse);
        cc = new CutComposer(prob, cuts, env, x_ij, logPath, nullptr);
        //
        cplex->setOut(env.getNullStream());
        if (logPath != "") {
            std::string _header("nItr,objV");
            for (CutBase *c: cc->cuts) {
                _header += "," + c->cut_name;
            }
            _header += ",note";
            char header[_header.size() + 1];
            std::strcpy(header, _header.c_str());
            createCSV(logPath, header);
        }
    }
    //
    ~RC() {
        delete cc;
    }
    rut::Solution* solve();
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
