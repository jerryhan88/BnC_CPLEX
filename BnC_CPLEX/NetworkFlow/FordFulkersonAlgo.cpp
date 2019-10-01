//
//  FordFulkersonAlgo.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 13/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "FordFulkersonAlgo.hpp"

void sendMaxFlow(double **rGraph, int numNodes, int source_id, int terminal_id, int parent[], double *max_flow) {
    int u, v;
    while (reachability_bfs(rGraph, numNodes, source_id, terminal_id, parent)) {
        // Find minimum residual capacity of the edges along the
        // path filled by BFS. Or we can say find the maximum flow
        // through the path found.
        double path_flow = DBL_MAX;
        for (v = terminal_id; v != source_id; v = parent[v]) {
            u = parent[v];
            if (rGraph[u][v] < path_flow) {
                path_flow = rGraph[u][v];
            }
        }
        // update residual capacities of the edges and reverse edges
        // along the path
        for (v = terminal_id; v != source_id; v = parent[v]) {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }

        // Add path flow to overall flow
        *max_flow += path_flow;
    }
}

std::set<std::set<int>> solve_multipleMaxFlow(double** x_ij, Problem* prob) {
    std::set<std::set<int>> outputSets;
    int V = (int) prob->N.size();
    double **rGraph = gen_graphInstance(x_ij, V);
    int parent[V];
    double max_flow;
    //
    bool isFirst = true;
    for (int i: prob->S) {
        if (i == prob->o || i == prob->d)
            continue;
        for (int j: prob->S) {
            if (i == j)
                continue;
            if (!isFirst) {
                for (int i = 0; i < V; i++) {
                    for (int j = 0; j < V; j++) {
                        rGraph[i][j] = x_ij[i][j];
                    }
                }
            } else {
                isFirst = false;
            }
            for (int i = 0; i < V; i++) {
                parent[i] = -1;
            }
            max_flow = 0.0;
            sendMaxFlow(rGraph, V, i, j, parent, &max_flow);
            bool visited[V];
            std::memset(visited, false, sizeof(visited));
            reachability_dfs(rGraph, V, i, visited);
            std::set<int> S1;
            for (int i = 0; i < V; i++) {
                if(visited[i]) {
                    S1.insert(i);
                }
            }
            outputSets.insert(S1);
        }
    }
    for (int i = 0; i < V; i++) {
        delete [] rGraph[i];
    }
    delete [] rGraph;
    //
    return outputSets;
}

