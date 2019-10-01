//
//  Problem.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 11/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Problem_hpp
#define Problem_hpp

#include <fstream>
#include <string>
#include <vector>
#include <set>

#include "nlohmann/json.hpp"


class Problem {
public:
    std::string problemName;
    int bv, bw, bu;
    //
    std::vector<int> K, PD, N;
    std::set<int> S, P, D;
    int o, d;
    int *h_k, *n_k, **c_ij;
    int *r_k, *v_k, *w_k;
    int *al_i, *be_i;
    int **t_ij;
    int M;
    //
    Problem () {};
    Problem(int numTasks,
                int *reward, int *volume, int *weight,
            int numNodes,
                int **distance, int **timeWindow,
            int bv, int bw, int bu);
    ~Problem();
    //
    void write_json(std::string ofpath);
    static Problem* read_json(std::string ifpath) {
        std::ifstream is(ifpath);
        nlohmann::json prob_json;
        is >> prob_json;
        is.close();
        //
        Problem *prob = new Problem();
        prob->problemName = prob_json["problemName"];
        for (int k: prob_json["K"]) {
            prob->K.push_back(k);
        }
        for (int i: prob_json["PD"]) {
            prob->PD.push_back(i);
        }
        for (int i: prob_json["N"]) {
            prob->N.push_back(i);
        }
        for (int i: prob_json["S"]) {
            prob->S.insert(i);
        }
        for (int i: prob_json["P"]) {
            prob->P.insert(i);
        }
        for (int i: prob_json["D"]) {
            prob->D.insert(i);
        }
        prob->o = prob_json["o"];
        prob->d = prob_json["d"];
        prob->r_k = new int[prob->K.size()];
        prob->v_k = new int[prob->K.size()];
        prob->w_k = new int[prob->K.size()];
        prob->h_k = new int[prob->K.size()];
        prob->n_k = new int[prob->K.size()];
        for (int k: prob->K) {
            prob->r_k[k] = prob_json["r_k"][k];
            prob->v_k[k] = prob_json["v_k"][k];
            prob->w_k[k] = prob_json["w_k"][k];
            prob->h_k[k] = prob_json["h_k"][k];
            prob->n_k[k] = prob_json["n_k"][k];
        }
        size_t numNodes = prob_json["t_ij"][0].size();
        prob->t_ij = new int*[numNodes];
        prob->c_ij = new int*[numNodes];
        for (int i = 0; i < numNodes; i++) {
            prob->t_ij[i] = new int[numNodes];
            prob->c_ij[i] = new int[numNodes];
        }
        prob->al_i = new int[numNodes];
        prob->be_i = new int[numNodes];
        for (int i: prob->N) {
            prob->al_i[i] = prob_json["al_i"][i];
            prob->be_i[i] = prob_json["be_i"][i];
            assert(0 <= prob->be_i[i]);
            for (int j: prob->N) {
                prob->t_ij[i][j] = prob_json["t_ij"][i][j];
                prob->c_ij[i][j] = prob_json["c_ij"][i][j];
            }
        }
        prob->M = prob_json["M"];
        prob->bv = prob_json["bv"];
        prob->bw = prob_json["bw"];
        prob->bu = prob_json["bu"];

        return prob;
    }
};

#endif /* Problem_hpp */
