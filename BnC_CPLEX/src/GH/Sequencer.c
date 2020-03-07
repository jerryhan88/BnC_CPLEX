//
//  Sequencer.c
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 31/10/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Sequencer.h"



void set_min_tt_seq(int* partialSeq, int partialSeqSize,
                     int* minSeq, int* seq,
                     int minSeqSize,
                     int seqBeginIndex4Search,
                     int n0, int n1,
                     double* al_i, double* be_i,
                     double** t_ij) {
    double min_tt = DBL_MAX;
    int i;
    for (i = 0; i < minSeqSize; i++) {
        minSeq[i] = -1;
    }
    if (minSeqSize - partialSeqSize == 1) {
        int s1;
        for (s1 = seqBeginIndex4Search; s1 < partialSeqSize; s1++) {
            for (i = 0; i < s1; i++) {
                seq[i] = partialSeq[i];
            }
            seq[s1] = n1;
            for (i = s1 + 1; i < minSeqSize; i++) {
                seq[i] = partialSeq[i - 1];
            }
            double tt = get_travelTime(seq, minSeqSize,
                                       al_i, be_i,
                                       t_ij);
            if (tt == -1.0) {
                continue;
            } else if (tt < min_tt) {
                if ( valid_TWs(seq, minSeqSize,
                                al_i, be_i,
                               t_ij) ) {
                    min_tt = tt;
                    for (i = 0; i < minSeqSize; i++) {
                        minSeq[i] = seq[i];
                    }
                }
            }
        }
    } else {
        int s0;
        for (s0 = 1; s0 < partialSeqSize; s0++) {
            int s1 = s0 > seqBeginIndex4Search ? s0 : seqBeginIndex4Search;
            for ( ; s1 < partialSeqSize; s1++) {
                for (i = 0; i < s0; i++) {
                    seq[i] = partialSeq[i];
                }
                seq[s0] = n0;
                for (i = s0 + 1; i < s1 + 1; i++) {
                    seq[i] = partialSeq[i - 1];
                }
                seq[s1 + 1] = n1;
                for (i = s1 + 2; i < minSeqSize; i++) {
                    seq[i] = partialSeq[i - 2];
                }
                double tt = get_travelTime(seq, minSeqSize,
                                            al_i, be_i,
                                            t_ij);
                if (tt == -1.0) {
                    continue;
                } else if (tt < min_tt) {
                    if ( valid_TWs(seq, minSeqSize,
                                    al_i, be_i,
                                   t_ij) ) {
                        min_tt = tt;
                        for (i = 0; i < minSeqSize; i++) {
                            minSeq[i] = seq[i];
                        }
                    }
                }
            }
        }
    }
}

bool valid_TWs(int* seq, int seqSize,
                double* al_i, double* be_i,
                double** t_ij) {
    int n0, n1;
    n0 = seq[0];
    double erest_deptTime, erest_arrvTime = -1.0;
    erest_deptTime = al_i[n0];
    int i;
    for (i = 1; i < seqSize; i++) {
        n1 = seq[i];
        erest_arrvTime = erest_deptTime + t_ij[n0][n1];
        if (be_i[n1] < erest_arrvTime) {
            return false;
        } else {
            erest_deptTime = erest_arrvTime > al_i[n1] ? erest_arrvTime : al_i[n1];
        }
        n0 = n1;
    }
    return true;
}

double get_travelTime(int* seq, int seqSize,
                      double* al_i, double* be_i,
                      double** t_ij) {
    int n0, n1;
    n0 = seq[0];
    double travelTime = 0.0;
    int i;
    for (i = 1; i < seqSize; i++) {
        n1 = seq[i];
        travelTime += t_ij[n0][n1];
        n0 = n1;
    }
    return travelTime;
}

double get_erest_arrvTime(int* seq, int seqSize,
                          double* al_i, double* be_i,
                          double** t_ij) {
    int n0, n1;
    n0 = seq[0];
    double erest_deptTime, erest_arrvTime = -1.0;
    erest_deptTime = al_i[n0];
    int i;
    for (i = 1; i < seqSize; i++) {
        n1 = seq[i];
        erest_arrvTime = erest_deptTime + t_ij[n0][n1];
        if (be_i[n1] < erest_arrvTime) {
            return -1.0;
        } else {
            erest_deptTime = erest_arrvTime > al_i[n1] ? erest_arrvTime : al_i[n1];
        }
        n0 = n1;
    }
    return erest_arrvTime;
}
