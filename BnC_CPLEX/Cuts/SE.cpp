//
//  SE.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 14/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//


#include "Base.hpp"

bool SE_cut::valid_subset(const std::set<int>& S1, Problem *prob) {
    if ( S1.size() <= 1 ||
        generatedSets.find(S1) != generatedSets.end() ) {
        return false;
    }
    std::set<int> setDiffRR;
    std::set_difference((*prob).S.begin(), (*prob).S.end(), S1.begin(), S1.end(), std::inserter(setDiffRR, setDiffRR.end()));
    if (setDiffRR.size() == 0) {
        std::set<int> vistedD;
        std::set_intersection((*prob).D.begin(), (*prob).D.end(), S1.begin(), S1.end(), std::inserter(vistedD, vistedD.end()));
        bool isAllPsVisited = true;
        for (int k: (*prob).K) {
            if (vistedD.find((*prob).n_k[k]) == vistedD.end()) {
                continue;
            }
            if (S1.find((*prob).h_k[k]) == S1.end()) {
                isAllPsVisited = false;
                break;
            }
        }
        return !isAllPsVisited;
    } else {
        return true;
    }
}

double get_LHS_SE0(const std::set<int>& S, double** x_ij, Problem* prob) {
    double lhs = 0.0;
    for (int i: S) {
        for (int j: S) {
            lhs += x_ij[i][j];
        }
    }
    lhs -= ((int) S.size() - 1);
    return lhs;
}

std::set<std::set<int>> solve_maxSE0(double** x_ij, Problem* prob) {
    std::set<std::set<int>> outputSets;
    run_TB4CG(x_ij, prob, get_LHS_SE0, outputSets);
    return outputSets;
}

SE_cut::SE_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd) : CutBase(ch_name, cutManagerType, isLocalCutAdd) {
    if (ch_name.substr(0, 1) == "e") {
        separationAlgo = solve_multipleMaxFlow;
    } else {
        assert(ch_name.substr(0, 1) == "h");
        separationAlgo = solve_maxSE0;
    }
}


std::set<std::set<int>> SE_cut::solve_separationProb(double** x_ij, CutComposer* cc) {
    Problem *prob = cc->prob;
    std::set<std::set<int>> validSets;
    std::set<std::set<int>> outputSets = separationAlgo(x_ij, prob);
    for (std::set<int> S1: outputSets) {
        if (valid_subset(S1, prob)) {
            double lhs = get_LHS_SE0(S1, x_ij, prob);
            if (lhs > 0) {
                validSets.insert(S1);
            }
        }
    }
    return validSets;
}

void SE_cut::set_LHS_Expr(IloExpr& lhs_expr, IloNumVar** x_ij, const std::set<int>& S1) {
    for (int i: S1) {
        for (int j: S1) {
            lhs_expr += x_ij[i][j];
        }
    }
    lhs_expr -= ((int) S1.size() - 1);
}

IloRangeArray SE_cut::get_cut_cnsts(double** x_ij, CutComposer* cc) {
    std::set<std::set<int>> validSets = solve_separationProb(x_ij, cc);
    //
    char buf[2048];
    IloRangeArray cnsts(cc->env);
    for (std::set<int> S1: validSets) {
        sprintf(buf, "%s(%d)", cut_name.c_str(), (int) generatedSets.size());
        generatedSets.insert(S1);
        IloExpr lhs_expr(cc->env);
        set_LHS_Expr(lhs_expr, cc->x_ij, S1);
        cnsts.add(lhs_expr <= 0);
        cnsts[cnsts.getSize() - 1].setName(buf);
        lhs_expr.end();
    }
    return cnsts;
}

void SE_cut::add_cnsts2Model(const std::set<std::set<int>>& validSets,  CutComposer* cc, const IloCplex::Callback::Context& context) {
    for (std::set<int> S1: validSets) {
        generatedSets.insert(S1);
        IloExpr lhs_expr(cc->env);
        set_LHS_Expr(lhs_expr, cc->x_ij, S1);
        addUserCutwCust(context, lhs_expr);
        lhs_expr.end();
        cc->numGenCuts += 1;
    }
}

void SE_cut::add_cut(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context) {
    std::set<std::set<int>> validSets = solve_separationProb(x_ij, cc);
    add_cnsts2Model(validSets, cc, context);
}

std::string SE_cut::add_cut_wLogging(double** x_ij, CutComposer* cc, const IloCplex::Callback::Context& context) {
    std::set<std::set<int>> validSets = solve_separationProb(x_ij, cc);
    add_cnsts2Model(validSets, cc, context);
    //
    std::string addedCuts;
    for (std::set<int> S1: validSets) {
        addedCuts += "(";
        for (int i: S1) {
            addedCuts += std::to_string(i) + "-";
        }
        addedCuts += ");";
    }
    return addedCuts;
}




/*

 ****************************
 Maybe used the following codes later!
 ****************************
 

 #define FROM_P1_PROB 45
 #define FROM_P2_PROB 90
 #define TASK_POP_RATIO 2
 #define TASK_GEN_RATIO 2
 #define SURVIVAL_RATIO 0.1
 
 
 void findOpt_candiSet(double** x_ij, Problem* prob,
                       std::set<std::set<int>>& validSets,
                       const std::set<std::set<int>>& generatedSets) {
     std::set<std::set<int>> outputSets = solve_multipleMaxFlow(x_ij, prob);
     for (std::set<int> S1: outputSets) {
         if (calc_LHS(x_ij, S1, prob, generatedSets) > 0) {
             validSets.insert(S1);
         }
     }
 }
 
 


class Individual {
public:
    int* chromosome;
    size_t c_size;
    std::set<int> S1;
    double fitness;
    Individual(int* chromosome, size_t c_size);
    ~Individual();
    //
    Individual* mate(Individual* parent2);
};

Individual::Individual(int* chromosome, size_t c_size) {
    this->chromosome = chromosome;
    this->c_size = c_size;
    for (int i = 0; i < this->c_size; i++) {
        if (this->chromosome[i] == 1) {
            S1.insert(i);
        }
    }
    fitness = -1.0;
}

Individual::~Individual() {
    delete [] chromosome;
}

Individual* Individual::mate(Individual* parent2) {
    int* child_chromosome = new int[c_size];
    for (int i = 0; i < c_size; i++) {
        int p1 = rand() % 100;
        if (p1 < FROM_P1_PROB) {
            child_chromosome[i] = chromosome[i];
        } else if (p1 < FROM_P2_PROB) {
            child_chromosome[i] = parent2->chromosome[i];
        } else {
            int p2 = rand() % 100;
            child_chromosome[i] = p2 < 50 ? 0 : 1;
        }
    }
    return new Individual(child_chromosome, c_size);
}

bool compareInd(Individual* ind1, Individual* ind2) {
    return ind1->fitness > ind2->fitness;
}

void findGA_candiSet(double** x_ij, Problem* prob,
                        std::set<std::set<int>>& validSets,
                     const std::set<std::set<int>>& generatedSets) {
    size_t population_size = (*prob).K.size() * TASK_POP_RATIO;
    std::vector<Individual*> population;
    for (int i = 0; i < population_size; i++) {
        int* chromosome = new int[(*prob).N.size()];
        for (int n = 0; n < (*prob).N.size(); n++) {
            int p = rand() % 100;
            chromosome[n] = p < 50 ? 0 : 1;
        }
        Individual* ind = new Individual(chromosome, (*prob).N.size());
        ind->fitness = calc_LHS(x_ij, ind->S1, prob, generatedSets);
        population.push_back(ind);
    }
    int generation = 0;
    while (generation < (*prob).K.size() * TASK_GEN_RATIO) {
        std::sort(population.begin(), population.end(), compareInd);
        std::vector<Individual*> new_generation;
        for (int i = 0; i < population.size(); i++) {
            if (population[i]->fitness <= 0)
                break;
            new_generation.push_back(population[i]);
        }
        size_t numSurvivors = new_generation.size();
        size_t numNewInds = population.size() - numSurvivors;
        for (int i = 0; i < numNewInds; i++) {
            int r = rand() % (population.size() / 2);
            Individual* parent1 = population[r];
            r = rand() % (population.size() / 2);
            Individual* parent2 = population[r];
            Individual* offspring = parent1->mate(parent2);
            offspring->fitness = calc_LHS(x_ij, offspring->S1, prob, generatedSets);
            new_generation.push_back(offspring);
        }
        
        for (int i = 0; i < numNewInds; i++) {
            delete population[numSurvivors + i];
        }
        population = new_generation;
        generation++;
    }
    for (Individual* ind: population) {
        if (ind->fitness > 0) {
            validSets.insert(ind->S1);
        }
        delete ind;
    }
}
*/
