//
//  TB.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 9/10/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef TB_hpp
#define TB_hpp

#include <set>
#include <deque>

#include "../Problem.hpp"
#include "../NetworkFlow/Graph.hpp"

void run_TB4CG(double** x_ij, Problem* prob,
               std::function<double (const std::set<int>& S, double** x_ij, Problem* prob)> get_LHS,
               std::set<std::set<int>>& outputSets);

#endif /* TB_hpp */
