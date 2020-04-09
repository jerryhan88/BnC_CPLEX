//
//  Solution.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 8/4/20.
//  Copyright © 2020 Chung-Kyun HAN. All rights reserved.
//

#include "../../include/ck_route/Other.hpp"


void rut::Solution::writeSolCSV(std::string solPathCSV) {
    char row[2048];
    std::fstream fout_csv;
    fout_csv.open(solPathCSV, std::ios::out);
    fout_csv << "objV,Gap,elaCpuTime,elaWallTime,notes" << "\n";
    sprintf(row, "%f,%f,%f,%f,%s",
            objV, gap, cpuT, wallT, note.c_str());
    fout_csv << row << "\n";
    fout_csv.close();
}

void rut::Solution::writeSolTXT(std::string solPathTXT) {
    std::fstream fout_txt;
    fout_txt.open(solPathTXT, std::ios::out);
    std::map<int, int> _route;
    for (int i: prob->N) {
        for (int j: prob->N) {
            if (x_ij[i][j] > 0.5) {
                _route[i] = j;
            }
        }
    }
    std::vector<int> seq;
    int i = prob->o;
    fout_txt << "Visiting sequence" << "\n";

    while (true) {
        fout_txt << i << "-";
        seq.push_back(i);
        i = _route[i];
        if (i == prob->d) {
            fout_txt << i << "\n\n";
            break;
        }
    }
    seq.push_back(i);
    fout_txt << "Arrival time" << "\n";
    for (int i: seq) {
        fout_txt << i << ": " << u_i[i] << "\n";
    }
    fout_txt.close();
}

