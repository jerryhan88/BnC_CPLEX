//
//  Sequencer.h
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 31/10/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Sequencer_h
#define Sequencer_h

#ifdef __cplusplus
extern "C" {
#endif


#include <stdio.h>
#include <float.h>
#include <stdbool.h>


void set_min_tt_seq(int* partialSeq, int partialSeqSize,
                     int* minSeq, int* seq,
                     int newSeqSize,
                     int seqBeginIndex4Search,
                     int n0, int n1,
                     double* al_i, double* be_i,
                     double** t_ij);

bool valid_TWs(int* seq, int seqSize,
                        double* al_i, double* be_i,
                        double** t_ij);


double get_travelTime(int* seq, int seqSize,
                      double* al_i, double* be_i,
                      double** t_ij);

double get_erest_arrvTime(int* seq, int seqSize,
                          double* al_i, double* be_i,
                          double** t_ij);

#ifdef __cplusplus
} //end extern "C"
#endif

#endif /* Sequencer_h */
