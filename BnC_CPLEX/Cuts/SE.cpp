//
//  SE.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 14/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//


#include "Base.hpp"
#include <float.h>

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


std::set<std::set<int>> solve_maxSE_constructCycle(double** x_ij, Problem* prob) {
    std::set<std::set<int>> outputSets;
    //
    int numNodes = (int) (*prob).N.size();
    std::vector<int> adjM_maxON;
    for (int i = 0; i < numNodes; i++) {
        int max_ON = -1;
        double max_OEV = -1.0;
        for (int j = 0; j < numNodes; j++) {
            if (i == j) {
                continue;
            }
            if (x_ij[i][j] > 0 && max_OEV < x_ij[i][j]){
                max_ON = j;
                max_OEV = x_ij[i][j];
            }
        }
        adjM_maxON.push_back(max_ON);
    }
    //
    bool visited[numNodes], isCycle;
    int n1;
    const int n2 = (*prob).d;
    std::vector<int> path;
    for (int n0: prob->S) {
        if (n0 == (*prob).o || n0 == (*prob).d) {
            continue;
        }
        std::memset(visited, false, sizeof(visited));
        isCycle = false;
        //
        path.clear();
        path.push_back(n0);
        while (n0 != n2) {
            visited[n0] = true;
            n1 = adjM_maxON[n0];
            if (n1 == -1) {
                break;
            }
            path.push_back(n1);
            if (visited[n1]) {
                isCycle = true;
                break;
            }
            n0 = n1;
        }
        if (isCycle) {
            std::set<int> S1(path.begin(), path.end());
            outputSets.insert(S1);
        }
    }
    return outputSets;
}

std::set<std::set<int>> solve_maxSE_DyBFS(double** x_ij, Problem* prob) {
    std::set<std::set<int>> outputSets;
    //
    int numNodes = (int) (*prob).N.size();
    
    int numPosFlow[numNodes];
    
//    std::cout << "\n";
    int** adjM;
    adjM = new int*[numNodes];
    for (int i = 0; i < numNodes; i++) {
        numPosFlow[i] = 0;
        adjM[i] = new int[numNodes];
        for (int j = 0; j < numNodes; j++) {
            if (x_ij[i][j] > 0.0 || x_ij[j][i] > 0.0) {
                adjM[i][numPosFlow[i]] = j;
                numPosFlow[i]++;
            }
//            std::cout << "x[" << i << "][" << j << "]: " << x_ij[i][j] << std::endl;
        }
    }
    double whole_LHS;
    std::set<int> S1;
//    std::cout << "\n";
    for (int n0: prob->S) {
        if (n0 == (*prob).o || n0 == (*prob).d) {
            continue;
        }
        whole_LHS = 0.0;
        S1.clear();
        S1.insert(n0);
        std::map<int, double> boundary;
        for (int i = 0; i < numPosFlow[n0]; i++) {
            int n1 = adjM[n0][i];
            double w = x_ij[n0][n1] + x_ij[n1][n0] - 1;
            boundary[n1] = w;
        }
//        std::cout << "\n";
        while (!boundary.empty()) {
            int max_n = -1;
            double max_w = -DBL_MAX;
            for (auto kv: boundary) {
                if (max_w < kv.second) {
                    max_n = kv.first;
                    max_w = kv.second;
                }
            }
            //
            assert(max_n != -1);
            S1.insert(max_n);
            whole_LHS += max_w;
//            std::cout << "(";
//            for (int i: S1) {
//                std::cout << i << ",";
//            }
//            std::cout << "): " << whole_LHS << std::endl;
            
            if (whole_LHS > 0.0) {
                outputSets.insert(S1);
                break;
            } else if (whole_LHS < -0.5) {
                break;
            }
            for (int i = 0; i < numPosFlow[max_n]; i++) {
                int n1 = adjM[max_n][i];
                if (S1.find(n1) == S1.end()) {
                    double w = 0.0;
                    boundary[n1] = w;
                }
            }
            boundary.erase(max_n);
            for (auto kv: boundary) {
                int n1 = kv.first;
                double w = -1.0;
                for (int n2: S1) {
                    w += x_ij[n2][n1] + x_ij[n1][n2];
                }
                boundary[n1] = w;
            }
        }
    }
    
    for (int i = 0; i < numNodes; i++) {
        delete [] adjM[i];
    }
    delete [] adjM;

    
    return outputSets;
}


SE_cut::SE_cut(std::string ch_name, IloCplex::CutManagement cutManagerType, IloBool isLocalCutAdd) : CutBase(ch_name, cutManagerType, isLocalCutAdd) {
    if (ch_name.substr(0, 1) == "e") {
        separationAlgo = solve_multipleMaxFlow;
    } else {
        assert(ch_name.substr(0, 1) == "h");
//        separationAlgo = solve_maxSE0;
        separationAlgo = solve_maxSE_constructCycle;
//        separationAlgo = solve_maxSE_DyBFS;
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
