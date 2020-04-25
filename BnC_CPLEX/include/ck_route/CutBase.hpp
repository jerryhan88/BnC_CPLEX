//
//  CutBase.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 1/3/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#ifndef CutBase_hpp
#define CutBase_hpp


#include <string>
#include <map>
#include <deque>
#include <mutex>
#include <cfloat>

#include <ilcplex/ilocplex.h>

#include "Other.hpp"

#include "ck_util/util.hpp"

typedef std::pair<int, int> edge;

class CutComposer;

int getNextNodeByFlow(int n0, CutComposer *cc, const IloCplex::Callback::Context &context);

int get_numP(std::set<int> S1, rut::Problem *prob);
double get_vS(const std::set<int> &S1, rut::Problem *prob);
double get_wS(const std::set<int> &S1, rut::Problem *prob);
double get_LHS_CA(const std::set<int> &S, double **x_ij, rut::Problem *prob);

int get_numFlowReturn(std::set<int> S1, rut::Problem *prob);

bool get_cyclicPath(int ori, int dest, bool *visited, int numNodes, std::vector<int> &path,
                    CutComposer *cc, const IloCplex::Callback::Context &context);

std::set<std::set<int>> get_maxCycles(const std::vector<int> &origins,
                                      CutComposer *cc, const IloCplex::Callback::Context &context);
std::set<std::set<edge>> get_infeasiblePaths(CutComposer *cc, const IloCplex::Callback::Context &context);


class CutBase {
public:
    std::string cut_name;
    IloCplex::CutManagement cutManagerType;
    IloBool isLocalCutAdd;
    int cutCounter;
    //
    CutBase(std::string cut_name){
        this->cut_name = cut_name;
    }
    CutBase(std::string cut_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd) : CutBase::CutBase(cut_name){
        this->cutManagerType = cutManagerType;
        this->isLocalCutAdd = isLocalCutAdd;
        this->cutCounter = 0;
    }
    //
    virtual IloRangeArray get_cut_cnsts(double **x_ij, CutComposer *cc) {
        throw "Should override get_cut_consts()";
    }
    virtual std::string add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context) {
        throw "Should override add_cut_wLogging()";
    }
    virtual void add_cut(CutComposer *cc, const IloCplex::Callback::Context &context) {
        throw "Should override add_cut()";
    }
    virtual IloRangeArray get_detectedCuts(CutComposer *cc) {
        throw "Should override get_detectedCuts()";
    }
    virtual void clear_detectedCuts() {
        throw "Should override clear_detectedCuts()";
    }
    //
    void addUserCutwCust(const IloCplex::Callback::Context &context, IloExpr lhs_expr);
};

class CutComposer : public IloCplex::Callback::Function {
public:
    rut::Problem *prob;
    std::vector<CutBase*> cuts;
    IloEnv env;
    IloNumVar **x_ij;
    double **_x_ij;
    std::string logPath;
    TimeTracker *tt;
    std::set<double> relaxedVals;
    int numGenCuts = 0;
    double time4Sep = 0.0;
    double time4FDV = 0.0;
    int num4Sep = 0;
    double bestRelVal = DBL_MAX;
    //
    CutComposer(rut::Problem *prob, std::vector<CutBase*> &cuts, IloEnv &env, IloNumVar **x_ij, std::string logPath, TimeTracker *tt);
    ~CutComposer() {
        for (CutBase* cut: cuts) {
            delete cut;
        }
        cuts.clear();
        for (int i: prob->N) {
            delete [] _x_ij[i];
        }
        delete [] _x_ij;
    };
    //
    virtual void invoke(const IloCplex::Callback::Context &context);
};


class SE_cut : public CutBase {
public:
    std::set<std::set<int>> generatedSets;
    SE_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd): CutBase(ch_name, cutManagerType, isLocalCutAdd) { }
    ~SE_cut() {
        generatedSets.clear();
    }
    //
    bool validate_subset(const std::set<int> &S1, rut::Problem *prob);
    bool validate_LHS(const std::set<int> &S1, CutComposer *cc, const IloCplex::Callback::Context &context);
    void add_cut(CutComposer *cc, const IloCplex::Callback::Context &context);
    void add_cnsts2Model(const std::set<std::set<int>> &validSets, CutComposer *cc, const IloCplex::Callback::Context &context);
    void clear_detectedCuts();
    void set_LHS_Expr(IloExpr &lhs_expr, IloNumVar **x_ij, const std::set<int> &S1);
    IloRangeArray get_cut_cnsts(double **x_ij, CutComposer *cc);
    //
    std::string add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context);
    IloRangeArray get_detectedCuts(CutComposer *cc);
private:
    std::set<std::set<int>> solve_separationProb(CutComposer *cc, const IloCplex::Callback::Context &context);
    
};

class CA_cut : public CutBase {
public:
    std::set<std::set<int>> generatedSets;
    CA_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd): CutBase(ch_name, cutManagerType, isLocalCutAdd) { }
    ~CA_cut() {
        generatedSets.clear();
    }
    //
    bool validate_subset(const std::set<int> &S1, rut::Problem *prob);
    bool validate_LHS(const std::set<int> &S1, rut::Problem *prob, CutComposer *cc, const IloCplex::Callback::Context &context);
    void add_cut(CutComposer *cc, const IloCplex::Callback::Context &context);
    void add_cnsts2Model(const std::set<std::set<int>> &validSets,  CutComposer *cc, const IloCplex::Callback::Context &context);
    void clear_detectedCuts();
    void set_LHS_Expr(IloExpr &lhs_expr, IloNumVar **x_ij, const std::set<int> &S1, rut::Problem *prob);
    IloRangeArray get_detectedCuts(CutComposer *cc);
    //
    std::string add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context);
    IloRangeArray get_cut_cnsts(double **x_ij, CutComposer *cc);
private:
    std::set<std::set<int>> solve_separationProb(CutComposer *cc, const IloCplex::Callback::Context &context);
};

// Routine Sequence cut!
class RS_cut : public CutBase {
public :
    std::set<std::set<int>> generatedSets;
    RS_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd): CutBase(ch_name, cutManagerType, isLocalCutAdd) { }
    ~RS_cut() {
        generatedSets.clear();
    }
    //
    bool validate_subset(const std::set<int> &S1, rut::Problem *prob);
    bool validate_LHS(const std::set<int> &S1, rut::Problem *prob, CutComposer *cc, const IloCplex::Callback::Context &context);
    void add_cut(CutComposer *cc, const IloCplex::Callback::Context &context);
    void add_cnsts2Model(const std::set<std::set<int>> &validSets, CutComposer *cc, const IloCplex::Callback::Context &context);
    void clear_detectedCuts();
    void set_LHS_Expr(IloExpr &lhs_expr, IloNumVar **x_ij, const std::set<int> &S1, rut::Problem *prob);
    IloRangeArray get_detectedCuts(CutComposer *cc);
    //
    std::string add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context);
    IloRangeArray get_cut_cnsts(double **x_ij, CutComposer *cc);
private:
    std::set<std::set<int>> solve_separationProb(CutComposer *cc, const IloCplex::Callback::Context &context);
};

// Infeasible Path cut!
class IP_cut : public CutBase {
public:
    std::set<std::set<edge>> generatedSets;
    IP_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd): CutBase(ch_name, cutManagerType, isLocalCutAdd) { }
    ~IP_cut() {
        generatedSets.clear();
    }
    //
    void add_cut(CutComposer *cc, const IloCplex::Callback::Context &context);
    void set_LHS_Expr(IloExpr &lhs_expr, IloNumVar **x_ij, const std::set<edge> &S1);
    IloRangeArray get_detectedCuts(CutComposer *cc);
    void clear_detectedCuts();
    //
    std::string add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context);
    IloRangeArray get_cut_cnsts(double **x_ij, CutComposer *cc);
private:
    bool validate_subset(const std::set<edge> &S1);
    std::set<std::set<edge>> solve_separationProb(CutComposer *cc, const IloCplex::Callback::Context &context);
    void add_cnsts2Model(const std::set<std::set<edge>> &validSets, CutComposer *cc, const IloCplex::Callback::Context &context);
};


class ALL_cut : public CutBase {
public:
    SE_cut* se_cut;
    CA_cut* ca_cut;
    RS_cut* rs_cut;
    IP_cut* ip_cut;
    ALL_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd): CutBase(ch_name, cutManagerType, isLocalCutAdd) {
        this->se_cut = new SE_cut("SE", cutManagerType, isLocalCutAdd);
        this->ca_cut = new CA_cut("CA", cutManagerType, isLocalCutAdd);
        this->rs_cut = new RS_cut("RS", cutManagerType, isLocalCutAdd);
        this->ip_cut = new IP_cut("IP", cutManagerType, isLocalCutAdd);
    }
    ~ALL_cut() {
        delete se_cut; delete ca_cut; delete rs_cut; delete ip_cut;
    }
    //
    void add_cut(CutComposer *cc, const IloCplex::Callback::Context &context);
    void clear_detectedCuts();
    IloRangeArray get_detectedCuts(CutComposer *cc);
    std::string add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context);
};

#endif /* CutBase_hpp */
