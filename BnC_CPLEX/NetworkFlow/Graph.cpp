//
//  Graph.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 13/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Graph.hpp"

double** gen_graphInstance(double **graph, int numNodes) {
    double **G;
    G = new double*[numNodes];
    for (int i = 0; i < numNodes; i++) {
        G[i] = new double[numNodes];
        for (int j = 0; j < numNodes; j++) {
            G[i][j] = graph[i][j];
        }
    }
    return G;
}

bool** gen_adjM(double **graph, int numNodes) {
    bool** adjM;
    adjM = new bool*[numNodes];
    for (int i = 0; i < numNodes; i++) {
        adjM[i] = new bool[numNodes];
        for (int j = 0; j < numNodes; j++) {
            if (graph[i][j] > 0.0) {
                adjM[i][j] = true;
            } else {
                adjM[i][j] = false;
            }
        }
    }
    return adjM;
}

bool reachability_bfs(double **rGraph, int V, int s, int t, int parent[]) {
    // Create a visited array and mark all vertices as not visited
    bool visited[V];
    std::memset(visited, 0, sizeof(visited));
    
    // Create a queue, enqueue source vertex and mark source vertex
    // as visited
    std::queue <int> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;
    
    // Standard BFS Loop
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        
        for (int v = 0; v < V; v++) {
            if (visited[v] == false && rGraph[u][v] > 0.0) {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }
    
    // If we reached sink in BFS starting from source, then return
    // true, else false
    return (visited[t] == true);
}

void reachability_dfs(double **graph, int V, int s, bool visited[]) {
    visited[s] = true;
    for (int i = 0; i < V; i++) {
        if (graph[s][i] && !visited[i])
            reachability_dfs(graph, V, i, visited);
    }
}
