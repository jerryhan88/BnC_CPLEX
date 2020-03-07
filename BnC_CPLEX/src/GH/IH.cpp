//
//  IH.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 31/10/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "IH.hpp"


InsertionHeuristic::InsertionHeuristic(Problem* prob, std::string logPath) {
    this->prob = prob;
    this->logPath = logPath;
    //
    x_ij = new double*[(*prob).N.size()];
    u_i = new double[(*prob).N.size()];
    for (int i: (*prob).N) {
        x_ij[i] = new double[(*prob).N.size()];
        for (int j: (*prob).N) {
            x_ij[i][j] = 0.0;
        }
        u_i[i] = 0.0;
    }
    
    //
    limDet = (*prob).bu;
    limVol = (*prob).bv;
    limWei = (*prob).bw;
    curReward = 0.0;
    curVol = 0.0;
    curWei = 0.0;
    //
    candiTasks.insert((*prob).K.begin(), (*prob).K.end());
    partialSeq.insert(partialSeq.end(), (*prob)._R.begin(), (*prob)._R.end());
    seqBeginIndex4Search = 1;
}

InsertionHeuristic::~InsertionHeuristic() {
    for (int i = 0; i < (*prob).N.size(); i++)
        delete [] x_ij[i];
    delete [] x_ij;
    delete [] u_i;
}

double get_travelTime(const std::vector<int>& seq, Problem* prob) {
    double* al_i = (*prob).al_i;
    double* be_i = (*prob).be_i;
    double** t_ij = (*prob).t_ij;
    //
    int n0, n1;
    n0 = seq[0];
    double erest_deptTime, erest_arrvTime = -1.0;
    erest_deptTime = al_i[n0];
    for (int i = 1; i < seq.size(); i++) {
        n1 = seq[i];
        erest_arrvTime = erest_deptTime + t_ij[n0][n1];
        if (be_i[n1] < erest_arrvTime) {
            return -1.0;
        } else {
            erest_deptTime = erest_arrvTime > al_i[n1] ? erest_arrvTime : al_i[n1];
        }
        n0 = n1;
    }
    assert(erest_arrvTime != -1.0);
    return erest_arrvTime;
}

tid_bestSeq InsertionHeuristic::findBestTaskSeq() {
    double best_reward = -1.0;
    std::vector<int> best_seq;
    int best_k = -1;
    double r, v, w;
    //
    int partialSeqSize, minSeqSize;
    partialSeqSize = (int) partialSeq.size();
    int* minSeq = new int[partialSeqSize + 2];
    int* seq = new int[partialSeqSize + 2];
    double* al_i = (*prob).al_i;
    double* be_i = (*prob).be_i;
    double** t_ij = (*prob).t_ij;
    
    for (int k: candiTasks) {
        r = (*prob).r_k[k];
        v = (*prob).v_k[k];
        w = (*prob).w_k[k];
        if (limVol < curVol + v) {
            continue;
        } else if (limWei < curWei + w) {
            continue;
        } else if (curReward + r <= best_reward) {
            continue;
        } else {
            if ( visitedWH.find( (*prob).h_k[k] ) != visitedWH.end() ) {
                minSeqSize = partialSeqSize + 1;
            } else {
                minSeqSize = partialSeqSize + 2;
            }
            int n0 = (*prob).h_k[k];
            int n1 = (*prob).n_k[k];
            set_min_tt_seq(&partialSeq[0], partialSeqSize,
                            minSeq, seq,
                            minSeqSize,
                            seqBeginIndex4Search,
                            n0, n1,
                            al_i, be_i,
                            t_ij);
            if (minSeq[0] == -1) {
                continue;
            }
            double tt = get_travelTime(minSeq, minSeqSize, al_i, be_i, t_ij);
            if (tt == -1.0) {
                continue;
            } else if (tt <= limDet) {
                if (best_reward < curReward + r) {
                    best_reward = curReward + r;
                    best_k = k;
                    best_seq.clear();
                    for (int i = 0; i < minSeqSize; i++) {
                        best_seq.push_back(minSeq[i]);
                    }
                }
            }
        }
    }
    //
    delete [] minSeq;
    delete [] seq;
    //
    return tid_bestSeq {best_k, best_seq};
}

void InsertionHeuristic::run() {
    while (candiTasks.size() != 0) {
        tid_bestSeq res = findBestTaskSeq();
        int best_k = res.first;
        if (best_k == -1) {
            break;
        }
        candiTasks.erase(best_k);
        visitedWH.insert((*prob).h_k[best_k]);
        partialSeq.clear();
        partialSeq.insert(partialSeq.end(), res.second.begin(), res.second.end());
        curVol += (*prob).v_k[best_k];
        curWei += (*prob).w_k[best_k];
        curReward += (*prob).r_k[best_k];
        int dp_index;
        for (dp_index = 0; dp_index < partialSeq.size(); dp_index++) {
            if (partialSeq[dp_index] == (*prob).n_k[best_k] ) {
                break;
            }
        }
        assert( dp_index != (int) partialSeq.size() );
        seqBeginIndex4Search = dp_index + 1;
    }
    //
    int n0, n1;
    n0 = partialSeq[0];
    double erest_deptTime, erest_arrvTime = (*prob).al_i[n0];
    u_i[n0] = erest_arrvTime;
    erest_deptTime = (*prob).al_i[n0];
    for (int i = 1; i < partialSeq.size(); i++) {
        n1 = partialSeq[i];
        erest_arrvTime = erest_deptTime + (*prob).t_ij[n0][n1];
        x_ij[n0][n1] = 1.0;
        erest_deptTime = erest_arrvTime > (*prob).al_i[n1] ? erest_arrvTime : (*prob).al_i[n1];
        n0 = n1;
        u_i[n0] = erest_arrvTime;
    }
    u_i[n0] = erest_arrvTime;
}

double InsertionHeuristic::get_objV() {
    double objV = 0.0;
    for (int k: (*prob).K) {
        for (int j: (*prob).N) {
            objV += (*prob).r_k[k] * x_ij[j][(*prob).n_k[k]];
        }
    }
    return objV;
}
