//
//  Graph.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 13/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Graph_hpp
#define Graph_hpp

#include <queue>
#include <cstring>

double** gen_graphInstance(double **graph, int numNodes);
//
bool reachability_bfs(double **rGraph, int V, int s, int t, int parent[]);
void reachability_dfs(double **graph, int V, int s, bool visited[]);


#endif /* Graph_hpp */
