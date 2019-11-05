//
//  Problem.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 11/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Problem.hpp"

Problem::~Problem() {
    std::vector<int*> ipV = {h_k, n_k};
    for(auto p: ipV) {
        delete p;
    }
    std::vector<double*> dpV = {r_k, v_k, w_k, al_i, be_i, v_i, w_i};
    for(auto p: dpV) {
        delete p;
    }
    
    for (int i = 0; i < N.size(); i++) {
        delete [] c_ij[i];
        delete [] t_ij[i];
    }
    delete [] c_ij;
    delete [] t_ij;
    //
    K.clear(); PD.clear(); S.clear(); N.clear();
    P.clear(); D.clear();
}


void Problem::write_json(std::string ofpath) {
    nlohmann::json prob_json;
    //
    prob_json["bv"] = bv;
    prob_json["bw"] = bw;
    prob_json["bu"] = bu;
    //
    prob_json["K"] = nlohmann::json (K);
    prob_json["PD"] = nlohmann::json (PD);
    prob_json["S"] = nlohmann::json (S);
    prob_json["N"] = nlohmann::json (N);
    prob_json["P"] = nlohmann::json (P);
    prob_json["D"] = nlohmann::json (D);
    prob_json["o"] = o;
    prob_json["d"] = d;
    std::vector<int> _r_k(K.size()), _v_k(K.size()), _w_k(K.size()), _h_k(K.size()), _n_k(K.size());
    for (int k: K) {
        _r_k[k] = r_k[k];
        _v_k[k] = v_k[k];
        _w_k[k] = w_k[k];
        _h_k[k] = h_k[k];
        _n_k[k] = n_k[k];
    }
    prob_json["r_k"] = nlohmann::json (_r_k);
    prob_json["v_k"] = nlohmann::json (_v_k);
    prob_json["w_k"] = nlohmann::json (_w_k);
    prob_json["h_k"] = nlohmann::json (_h_k);
    prob_json["n_k"] = nlohmann::json (_n_k);
    //
    std::vector<double> _al_i(N.size()), _be_i(N.size());
    std::vector<std::vector<double>> _t_ij(N.size());
    std::vector<std::vector<int>> _c_ij(N.size());
    for (int i: N) {
        _al_i[i] = al_i[i];
        _be_i[i] = be_i[i];
        std::vector<double> _t_i(N.size());
        std::vector<int> _c_i(N.size());
        for (int j: N) {
            _t_i[j] = t_ij[i][j];
            _c_i[j] = c_ij[i][j];
        }
        _t_ij[i] = _t_i;
        _c_ij[i] = _c_i;
    }
    prob_json["al_i"] = nlohmann::json (_al_i);
    prob_json["be_i"] = nlohmann::json (_be_i);
    prob_json["t_ij"] = nlohmann::json (_t_ij);
    prob_json["c_ij"] = nlohmann::json (_c_ij);
    //
    prob_json["M"] = M;

    std::ofstream os(ofpath);
    os << prob_json <<std::endl;

}
