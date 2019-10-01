//
//  Etc.hpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 11/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#ifndef Etc_hpp
#define Etc_hpp

#include <fstream>
#include <ctime>
#include <chrono>
#include <string>

void createCSV(std::string fpath, char *header);
void appendRow(std::string fpath, char *row);

class TimeTracker {
public:
    std::clock_t c_start;
    std::chrono::high_resolution_clock::time_point w_start;
    //
    TimeTracker() {
        c_start = std::clock();
        w_start = std::chrono::high_resolution_clock::now();
    }
    ~TimeTracker() {}
    //
    std::string get_curTime();
    double get_elipsedTimeCPU();
    double get_elipsedTimeWall();
};

class FilePathOrganizer {
public:
    std::string logPath, solPathCSV, solPathTXT, lpPath, ilpPath;
    //
    FilePathOrganizer(const std::string &appr_dpath, const std::string &postfix);
    ~FilePathOrganizer() {}
};



#endif /* Etc_hpp */
