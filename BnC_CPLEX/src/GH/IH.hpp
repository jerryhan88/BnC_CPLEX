//
//  IH.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 31/10/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef IH_hpp
#define IH_hpp

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <cfloat>

#include "../../include/ck_route/Problem.hpp"
#include "Sequencer.h"
//extern "C" {
//    #include "Sequencer.h"
//}

typedef std::pair< int, std::vector<int> > tid_bestSeq;

using namespace rut;

class InsertionHeuristic {
public:
    Problem* prob;
    std::string logPath;
    double** x_ij;
    double* u_i;
    double limDet, limVol, limWei;
    double curReward, curVol, curWei;
    std::set<int> candiTasks, visitedWH;
    std::vector<int> partialSeq;
    int seqBeginIndex4Search;
    //
    InsertionHeuristic(Problem* prob, std::string logPath);
    ~InsertionHeuristic();
    void run();
    tid_bestSeq findBestTaskSeq();
    double get_objV();
};


#endif /* IH_hpp */
