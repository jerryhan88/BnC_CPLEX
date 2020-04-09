//
//  IH.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 31/10/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ck_route/Router.hpp"

rut::Solution* rgh::InsertionHeuristic::solve() {
    run();
    //
    rut::Solution *sol = new rut::Solution(prob);
    sol->objV = get_objV();
    sol->gap = -1.0;
    sol->cpuT = tt->get_elapsedTimeCPU();
    sol->wallT = tt->get_elapsedTimeWall();
    
    for (int i: prob->N) {
        for (int j: prob->N) {
            sol->x_ij[i][j] = x_ij[i][j];
        }
        sol->u_i[i] = u_i[i];
    }
    return sol;
}

void rgh::InsertionHeuristic::run() {
    while (candiTasks.size() != 0) {
        tid_bestSeq res = findBestTaskSeq();
        int best_k = res.first;
        if (best_k == -1) {
            break;
        }
        candiTasks.erase(best_k);
        visitedWH.insert(prob->h_k[best_k]);
        partialSeq.clear();
        partialSeq.insert(partialSeq.end(), res.second.begin(), res.second.end());
        curVol += prob->v_k[best_k];
        curWei += prob->w_k[best_k];
        curReward += prob->r_k[best_k];
        int dp_index;
        for (dp_index = 0; dp_index < partialSeq.size(); dp_index++) {
            if (partialSeq[dp_index] == prob->n_k[best_k] ) {
                break;
            }
        }
        assert( dp_index != (int) partialSeq.size() );
        seqBeginIndex4Search = dp_index + 1;
    }
    //
    int n0, n1;
    n0 = partialSeq[0];
    double erest_deptTime, erest_arrvTime = prob->al_i[n0];
    u_i[n0] = erest_arrvTime;
    erest_deptTime = prob->al_i[n0];
    for (int i = 1; i < partialSeq.size(); i++) {
        n1 = partialSeq[i];
        erest_arrvTime = erest_deptTime + prob->t_ij[n0][n1];
        x_ij[n0][n1] = 1.0;
        erest_deptTime = erest_arrvTime > prob->al_i[n1] ? erest_arrvTime : prob->al_i[n1];
        n0 = n1;
        u_i[n0] = erest_arrvTime;
    }
    u_i[n0] = erest_arrvTime;
}

void rgh::InsertionHeuristic::initIH() {
    for (int i: prob->N) {
        for (int j: prob->N) {
            x_ij[i][j] = 0.0;
        }
        u_i[i] = 0.0;
    }
    curReward = 0.0;
    curVol = 0.0;
    curWei = 0.0;
    //
    candiTasks.clear();
    visitedWH.clear();
    partialSeq.clear();
    candiTasks.insert(prob->K.begin(), prob->K.end());
    partialSeq.insert(partialSeq.end(), prob->_R.begin(), prob->_R.end());
    seqBeginIndex4Search = 1;
}

tid_bestSeq rgh::InsertionHeuristic::findBestTaskSeq() {
    double best_reward = -1.0;
    std::vector<int> best_seq;
    int best_k = -1;
    double r, v, w;
    //
    int partialSeqSize, minSeqSize;
    partialSeqSize = (int) partialSeq.size();
    int *minSeq = new int[partialSeqSize + 2];
    int *seq = new int[partialSeqSize + 2];
    double *al_i = prob->al_i;
    double *be_i = prob->be_i;
    double **t_ij = prob->t_ij;
    
    for (int k: candiTasks) {
        r = prob->r_k[k];
        v = prob->v_k[k];
        w = prob->w_k[k];
        if (limVol < curVol + v) {
            continue;
        } else if (limWei < curWei + w) {
            continue;
        } else if (curReward + r <= best_reward) {
            continue;
        } else {
            if ( visitedWH.find( prob->h_k[k] ) != visitedWH.end() ) {
                minSeqSize = partialSeqSize + 1;
            } else {
                minSeqSize = partialSeqSize + 2;
            }
            int n0 = prob->h_k[k];
            int n1 = prob->n_k[k];
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
            if (tt <= limDet) {
                if ( valid_TWs(minSeq, minSeqSize,
                                al_i, be_i,
                               t_ij) ) {
                    if (best_reward < curReward + r) {
                        best_reward = curReward + r;
                        best_k = k;
                        best_seq.clear();
                        for (int i = 0; i < minSeqSize; i++) {
                            best_seq.push_back(minSeq[i]);
                        }
                    }
                } else {
                    continue;
                }
            } else {
                continue;
            }
        }
    }
    //
    delete [] minSeq;
    delete [] seq;
    //
    return tid_bestSeq {best_k, best_seq};
}

double rgh::InsertionHeuristic::get_objV() {
    double objV = 0.0;
    for (int k: prob->K) {
        for (int j: prob->N) {
            objV += prob->r_k[k] * x_ij[j][prob->n_k[k]];
        }
    }
    return objV;
}
