//
//  TB.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 9/10/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "TB.hpp"


#define TL_SIZE_RATIO 0.2
#define TB_ITER_RATIO 0.75

void set_adjN_P(const std::set<int>& S1, bool **adjM, size_t numNodes,
                std::set<int>& adjN_P) {
    adjN_P.clear();
    for (int i = 0; i < numNodes; i++) {
        if (S1.find(i) != S1.end()) {
            continue;
        }
        for (int j: S1) {
            if (adjM[i][j] || adjM[j][i]) {
                adjN_P.insert(i);
            }
        }
    }
}

void set_adjN_M(const std::set<int>& S1, bool **adjM, size_t numNodes,
             std::set<int>& adjN_M) {
    adjN_M.clear();
    for (int i: S1) {
        for (int j = 0; j < numNodes; j++) {
            if (S1.find(j) != S1.end()) {
                continue;
            }
            if (adjM[i][j] || adjM[j][i]) {
                adjN_M.insert(i);
            }
        }
    }
}

void run_TB4CG(double** x_ij, Problem* prob,
               std::function<double (const std::set<int>& S, double** x_ij, Problem* prob)> get_LHS,
               std::set<std::set<int>>& outputSets) {
    int numNodes = (int) prob->N.size();
    bool** adjM = gen_adjM(x_ij, numNodes);
    int tabu_list_size = (int) numNodes * TL_SIZE_RATIO;
    int ts_num_iteration = (int) numNodes * TB_ITER_RATIO;
    std::set<int> S1, S2;
    std::set<int> adjN_P, adjN_M;
    for (int i: prob->S) {
        S1.clear();
        S1.insert(i);
        //
        std::deque<int> tabu_list;
        int num_iter = 0;
        while (num_iter < ts_num_iteration) {
            double maxLHS = -1.0;
            int max_j = -1;
            set_adjN_P(S1, adjM, numNodes, adjN_P);
            for (int j: adjN_P) {
                std::deque<int>::iterator it = tabu_list.begin();
                bool isInTB = false;
                while(it != tabu_list.end()) {
                    if (j == *it) {
                        isInTB = true;
                        break;
                    }
                    it++;
                }
                if (isInTB) {
                    continue;
                }
                S2.clear();
                S2.insert(S1.begin(), S1.end());
                S2.insert(j);
                double lhs = get_LHS(S2, x_ij, prob);
                if (maxLHS < lhs) {
                    maxLHS = lhs;
                    max_j = j;
                }
            }
            set_adjN_M(S1, adjM, numNodes, adjN_M);
            for (int j: adjN_M) {
                if (i == j) {
                    continue;
                }
                S2.clear();
                S2.insert(S1.begin(), S1.end());
                S2.erase(j);
                if (S2.size() == 0) {
                    continue;
                }
                double lhs = get_LHS(S2, x_ij, prob);
                if (maxLHS < lhs) {
                    maxLHS = lhs;
                    max_j = j;
                }
            }
            if (max_j == -1) {
                break;
            }
            if (adjN_P.find(max_j) != adjN_P.end()) {
                S1.insert(max_j);
            } else {
                S1.erase(max_j);
                tabu_list.push_back(max_j);
                if ( tabu_list_size < tabu_list.size() ) {
                    tabu_list.pop_front();
                }
            }
            if (outputSets.find(S1) != outputSets.end()) {
                break;
            } else {
                outputSets.insert(S1);
            }
            num_iter++;
        }
    }
    for (int i = 0; i < numNodes; i++) {
        delete [] adjM[i];
    }
    delete [] adjM;
}


