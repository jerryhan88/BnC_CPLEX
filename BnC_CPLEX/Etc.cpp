//
//  Etc.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 11/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include "Etc.hpp"


void createCSV(std::string fpath, char *header) {
    if (fpath == "") {
        return;
    }
    std::fstream fout;
    fout.open(fpath, std::ios::out);
    fout << header << "\n";
    fout.close();
}

void appendRow(std::string fpath, char *row) {
    if (fpath == "") {
        return;
    }
    std::fstream fout;
    fout.open(fpath, std::ios::out | std::ios::app);
    fout << row << "\n";
    fout.close();
}

std::string TimeTracker::get_curTime() {
    time_t now = time(0);
    char* dt = ctime(&now);
    return std::string(dt);
}

double TimeTracker::get_elipsedTimeCPU() {
    std::clock_t c_end = std::clock();
    return (c_end-c_start) / (double) CLOCKS_PER_SEC;
}

double TimeTracker::get_elipsedTimeWall() {
    std::chrono::high_resolution_clock::time_point w_end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(w_end-w_start).count() / 1000.0;
}

FilePathOrganizer::FilePathOrganizer(const std::string &appr_dpath, const std::string &postfix) {
    std::string appr_name = appr_dpath.substr(appr_dpath.find_last_of("/") + 1, appr_dpath.size());
    std::string prefix(appr_dpath + "/" + appr_name + "_" + postfix);
    //
    this->logPath = prefix + ".log";
    this->solPathCSV = prefix + ".csv";
    this->solPathTXT = prefix + ".txt";
    this->lpPath = prefix + ".lp";
    this->ilpPath = prefix + ".ilp";
}
