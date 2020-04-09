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
    virtual void invoke (const IloCplex::Callback::Context &context);
};


class SE_cut : public CutBase {
public:
    std::set<std::set<int>> generatedSets;
    std::set<std::set<int>> addedSets2MM;
    std::function<std::set<std::set<int>>(CutComposer *cc, const IloCplex::Callback::Context &context)> separationAlgo;
    SE_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd);
    //
    IloRangeArray get_cut_cnsts(double **x_ij, CutComposer *cc);
    void add_cut(CutComposer *cc, const IloCplex::Callback::Context &context);
    std::string add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context);
    IloRangeArray get_detectedCuts(CutComposer *cc);
private:
    bool valid_subset(const std::set<int> &S1, rut::Problem *prob);
    std::set<std::set<int>> solve_separationProb(CutComposer *cc, const IloCplex::Callback::Context &context);
    void set_LHS_Expr(IloExpr &lhs_expr, IloNumVar **x_ij, const std::set<int> &S1);
    void add_cnsts2Model(const std::set<std::set<int>> &validSets, CutComposer *cc, const IloCplex::Callback::Context &context);
};

class CA_cut : public CutBase {
public:
    std::set<std::set<int>> generatedSets;
    std::set<std::set<int>> addedSets2MM;
    std::function<std::set<std::set<int>>(CutComposer *cc, const IloCplex::Callback::Context &context)> separationAlgo;
    CA_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd);
    //
    IloRangeArray get_cut_cnsts(double **x_ij, CutComposer *cc);
    void add_cut(CutComposer *cc, const IloCplex::Callback::Context &context);
    std::string add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context);
    IloRangeArray get_detectedCuts(CutComposer *cc);
private:
    bool valid_subset(const std::set<int> &S1, rut::Problem *prob);
    std::set<std::set<int>> solve_separationProb(CutComposer *cc, const IloCplex::Callback::Context &context);
    void set_LHS_Expr(IloExpr &lhs_expr, IloNumVar **x_ij, const std::set<int> &S1, rut::Problem *prob);
    void add_cnsts2Model(const std::set<std::set<int>> &validSets,  CutComposer *cc, const IloCplex::Callback::Context &context);
};

// Routine Sequence cut!
class RS_cut : public CutBase {
public :
    std::set<std::set<int>> generatedSets;
    std::set<std::set<int>> addedSets2MM;
    std::function<std::set<std::set<int>>(CutComposer *cc, const IloCplex::Callback::Context &context)> separationAlgo;
    RS_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd);
    IloRangeArray get_cut_cnsts(double **x_ij, CutComposer *cc);
    void add_cut(CutComposer *cc, const IloCplex::Callback::Context &context);
    std::string add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context);
    IloRangeArray get_detectedCuts(CutComposer *cc);
private:
    bool valid_subset(const std::set<int> &S1, rut::Problem *prob);
    std::set<std::set<int>> solve_separationProb(CutComposer *cc, const IloCplex::Callback::Context &context);
    void set_LHS_Expr(IloExpr &lhs_expr, IloNumVar **x_ij, const std::set<int> &S1, rut::Problem *prob);
    void add_cnsts2Model(const std::set<std::set<int>> &validSets, CutComposer *cc, const IloCplex::Callback::Context &context);
};

// Infeasible Path cut!
class IP_cut : public CutBase {
public:
    std::set<std::set<edge>> generatedSets;
    std::set<std::set<edge>> addedSets2MM;
    std::function<std::set<std::set<edge>>(CutComposer *cc, const IloCplex::Callback::Context &context)> separationAlgo;
    IP_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd);
    //
    IloRangeArray get_cut_cnsts(double **x_ij, CutComposer *cc);
    void add_cut(CutComposer *cc, const IloCplex::Callback::Context &context);
    std::string add_cut_wLogging(CutComposer *cc, const IloCplex::Callback::Context &context);
    IloRangeArray get_detectedCuts(CutComposer *cc);
private:
    bool valid_subset(const std::set<edge> &S1);
    std::set<std::set<edge>> solve_separationProb(CutComposer *cc, const IloCplex::Callback::Context &context);
    void set_LHS_Expr(IloExpr &lhs_expr, IloNumVar **x_ij, const std::set<edge> &S1);
    void add_cnsts2Model(const std::set<std::set<edge>> &validSets, CutComposer *cc, const IloCplex::Callback::Context &context);
};

#endif /* CutBase_hpp */
