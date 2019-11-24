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
