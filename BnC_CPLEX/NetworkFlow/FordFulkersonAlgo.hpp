//
//  FordFulkersonAlgo.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 13/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef FordFulkersonAlgo_hpp
#define FordFulkersonAlgo_hpp

#include <set>
#include <cfloat>

#include "Graph.hpp"
#include "../Problem.hpp"


void sendMaxFlow(double **rGraph, int numNodes, int source_id, int terminal_id, int parent[], double *max_flow);

std::set<std::set<int>> solve_multipleMaxFlow(double** x_ij, Problem* prob);

#endif /* FordFulkersonAlgo_hpp */
