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
    double bv, bw, bu;
    //
    std::vector<int> K, PD, N, _S;
    std::vector<std::set<int>> LS;
    std::set<int> P, D, S;
    int o, d;
    int *h_k, *n_k, **c_ij;
    double *r_k, *v_k, *w_k;
    double *al_i, *be_i;
    double **t_ij;
    double M;
    double *v_i, *w_i;
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
            prob->_S.push_back(i);
        }
        for (int i = 0; i < prob->_S.size() - 1; i++) {
            std::set<int> lS;
            lS.insert(prob->_S.begin() + i + 2, prob->_S.end());
            prob->LS.push_back(lS);
        }
        for (int i: prob_json["P"]) {
            prob->P.insert(i);
        }
        std::vector<int> _D;
        for (int i: prob_json["D"]) {
            prob->D.insert(i);
            _D.push_back(i);
        }
        prob->o = prob_json["o"];
        prob->d = prob_json["d"];
        prob->r_k = new double[prob->K.size()];
        prob->v_k = new double[prob->K.size()];
        prob->w_k = new double[prob->K.size()];
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
        prob->t_ij = new double*[numNodes];
        prob->c_ij = new int*[numNodes];
        for (int i = 0; i < numNodes; i++) {
            prob->t_ij[i] = new double[numNodes];
            prob->c_ij[i] = new int[numNodes];
        }
        prob->al_i = new double[numNodes];
        prob->be_i = new double[numNodes];
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
        //
        prob->v_i = new double[numNodes];
        prob->w_i = new double[numNodes];
        for (int n0: prob->S) {
            prob->v_i[n0] = 0.0;
            prob->w_i[n0] = 0.0;
        }
        for (int n0: prob->P) {
            prob->v_i[n0] = 0.0;
            prob->w_i[n0] = 0.0;
        }
        
        for (int k = 0; k < _D.size(); k++) {
            prob->v_i[_D[k]] = prob->v_k[k];
            prob->w_i[_D[k]] = prob->w_k[k];
        }
        return prob;
    }
};

#endif /* Problem_hpp */
