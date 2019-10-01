//
//  Base.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 13/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Base_hpp
#define Base_hpp

#include <string>
#include <map>
#include <mutex>

#include <ilcplex/ilocplex.h>

#include "../Problem.hpp"
#include "../Etc.hpp"



class CutComposer;

class CutBase {
public:
    std::string cut_name;
    //
    CutBase(std::string cut_name){
        this->cut_name = cut_name;
    }
    //
    virtual IloRangeArray get_cut_cnsts(double** x_ij, CutComposer* cc) {
        throw "Should override get_cut_consts()";
    }
    virtual std::string add_cut_wLogging(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context) {
        throw "Should override add_cut_wLogging()";
    }
    virtual void add_cut(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context) {
        throw "Should override add_cut()";
    }
};

class CutComposer : public IloCplex::Callback::Function {
public:
    Problem* prob;
    std::vector<CutBase*> cuts;
    IloEnv env;
    IloNumVar** x_ij;
    std::string logPath;
    //
    CutComposer(Problem* prob, std::vector<CutBase*>& cuts, IloEnv& env, IloNumVar** x_ij, std::string logPath);
    ~CutComposer() {
        for (CutBase* cut: cuts) {
            delete cut;
        }
        cuts.clear();
    };
    //
    virtual void invoke (const IloCplex::Callback::Context &context);
    //
private:
    double** get_x_ij(const IloCplex::Callback::Context &context);
    void add_cuts();
};


typedef std::pair<int, int> edge;



class SE0_cut : public CutBase {
public:
    std::set<std::set<int>> generatedSets;
//    std::function<void(double** x_ij, Problem* prob, std::set<std::set<int>>& validSets, const std::set<std::set<int>>& generatedSets)> find_candiSet;
    SE0_cut(std::string ch_name);
    //
    IloRangeArray get_cut_cnsts(double** x_ij, CutComposer* cc);
    void add_cut(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context);
    std::string add_cut_wLogging(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context);
private:
    bool valid_subset(const std::set<int>& S1, Problem *prob);
    std::set<std::set<int>> solve_separationProb(double** x_ij, CutComposer* cc);
    void set_LHS_Expr(IloExpr& lhs_expr, IloNumVar** x_ij, const std::set<int>& S1);
    void add_cnsts2Model(const std::set<std::set<int>>& validSets, CutComposer* cc, const IloCplex::Callback::Context& context);
};

class Capacity_cut : public CutBase {
public:
    std::set<std::set<int>> generatedSets;
    Capacity_cut(std::string ch_name) : CutBase(ch_name) {}
    //
    IloRangeArray get_cut_cnsts(double** x_ij, CutComposer* cc);
    void add_cut(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context);
    std::string add_cut_wLogging(double **x_ij, CutComposer *cc, const IloCplex::Callback::Context &context);
private:
    bool valid_subset(const std::set<int>& S1, Problem *prob);
    std::set<std::set<int>> solve_separationProb(double** x_ij, CutComposer* cc);
    void set_LHS_Expr(IloExpr& lhs_expr, IloNumVar** x_ij, const std::set<int>& S1, Problem* prob);
    void add_cnsts2Model(const std::set<std::set<int>>& validSets,  CutComposer* cc, const IloCplex::Callback::Context& context);
};

class Cover_cut : public CutBase {
public:
    std::set<std::set<edge>> generatedSets;
    Cover_cut(std::string ch_name) : CutBase(ch_name) {}
    //
    IloRangeArray get_cut_cnsts(double** x_ij, CutComposer* cc);
    void add_cut(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context);
    std::string add_cut_wLogging(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context);
private:
    bool valid_subset(const std::set<edge>& S1);
    std::set<edge> solve_separationProb(double** x_ij, CutComposer* cc);
    void set_LHS_Expr(IloExpr& lhs_expr, IloNumVar** x_ij, const std::set<edge>& S1);
    void add_cnsts2Model(const std::set<edge>& S1, CutComposer* cc, const IloCplex::Callback::Context& context);
};

#endif /* Base_hpp */
