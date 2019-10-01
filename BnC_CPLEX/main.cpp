//
//  main.cpp
//  BnC_CPLEX
//
//  Created by Chung-Kyun HAN on 11/9/19.
//  Copyright © 2019 Chung-Kyun HAN. All rights reserved.
//

#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <dirent.h>

#include "Problem.hpp"
#include "Etc.hpp"
#include "MathematicalModel/BaseMM.hpp"


#define DEFAULT_BUFFER_SIZE 2048

bool hasOption(std::vector<std::string> &arguments, std::string option) {
    bool hasValue = false;
    for (std::string str: arguments) {
        if (str == option) {
            hasValue = true;
        }
    }
    return hasValue;
}

std::string valueOf(std::vector<std::string> &arguments, std::string option) {
    std::string value = "";
    for (int i = 0; i < arguments.size(); i++) {
        if (arguments[i] == option) {
            value = arguments[i + 1];
        }
    }
    return value;
}

std::vector<std::string> parseWithDelimiter(std::string str, std::string delimiter) {
    std::vector<std::string> tokens;
    size_t pos = 0;
    std::string token;
    while ((pos = str.find(delimiter)) != std::string::npos) {
        token = str.substr(0, pos);
        tokens.push_back(token);
        str.erase(0, pos + delimiter.length());
    }
    tokens.push_back(str);
    return tokens;
    
}

std::vector<std::string> read_directory(const std::string &d_path, const std::string &extension) {
    std::vector<std::string> fileNames;
    DIR* dirp = opendir(d_path.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        std::string fn(dp->d_name);
        std::size_t found = fn.find(extension);
        if (found!=std::string::npos)
            fileNames.push_back(fn);
    }
    closedir(dirp);
    std::sort(fileNames.begin(), fileNames.end());
    return fileNames;
}

void write_solution(Problem *prob, FilePathOrganizer &fpo, TimeTracker &tt, BaseMM *mm) {
    char row[DEFAULT_BUFFER_SIZE];
    std::fstream fout_csv;
    fout_csv.open(fpo.solPathCSV, std::ios::out);
    fout_csv << "objV,eliCpuTime,eliWallTime" << "\n";
    sprintf(row, "%f,%f,%f",
            mm->cplex->getObjValue(),
            tt.get_elipsedTimeCPU(),
            tt.get_elipsedTimeWall());
    fout_csv << row << "\n";
    fout_csv.close();
}

void write_solution(Problem *prob, FilePathOrganizer &fpo, TimeTracker &tt,
                    BaseMM *mm,
                    double gap) {
    char row[DEFAULT_BUFFER_SIZE];
    std::fstream fout_csv;
    fout_csv.open(fpo.solPathCSV, std::ios::out);
    fout_csv << "objV,Gap,eliCpuTime,eliWallTime" << "\n";
    sprintf(row, "%f,%f,%f,%f",
            mm->cplex->getObjValue(),
            gap,
            tt.get_elipsedTimeCPU(),
            tt.get_elipsedTimeWall());
    fout_csv << row << "\n";
    fout_csv.close();
    //
    double **x_ij = new double *[(*prob).N.size()];
    double *u_i = new double[(*prob).N.size()];
    for (int i: (*prob).N) {
        x_ij[i] = new double[(*prob).N.size()];
    }
    mm->get_x_ij(x_ij);
    mm->get_u_i(u_i);
    std::map<int, int> _route;
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            if (x_ij[i][j] > 0.5) {
                _route[i] = j;
            }
        }
    }
    
    std::fstream fout_txt;
    fout_txt.open(fpo.solPathTXT, std::ios::out);
    
    
//    fout_txt <<"\n";
//
//    for (int k: (*prob).K) {
//        double r = 0.0;
//        for (int j: (*prob).N) {
//            r += (*prob).r_k[k] * x_ij[j][(*prob).n_k[k]];
//        }
//        if (r > 0.5) {
//            fout_txt << k << "\n";
//        }
//    }
//    fout_txt <<"\n";
//
//    for (int i: (*prob).N) {
//        for (int j: (*prob).N) {
//            if (x_ij[i][j] > 0.5) {
//                fout_txt << "x[" << i << "][" << j << "]:" << x_ij[i][j];
//                fout_txt << "(" << (*prob).t_ij[i][j] << ")" << "\n";
//            }
//        }
//    }
//    fout_txt <<"\n";
    
    std::vector<int> seq;
    int i = (*prob).o;
    fout_txt << "Visiting sequence" << "\n";

    while (true) {
        fout_txt << i << "-";
        seq.push_back(i);
        i = _route[i];
        if (i == (*prob).d) {
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
    //

    for (int i: (*prob).N) {
        delete [] x_ij[i];
    }
    delete [] x_ij;
    delete [] u_i;
}


int main(int argc, const char * argv[]) {
    std::vector<std::string> arguments;
    for (int i = 0; i < argc; i++) {
        arguments.push_back(argv[i]);
    }
    std::string basicFlags[] = {"-i", "-a"};
    int numBasicFlags = sizeof(basicFlags) / sizeof(std::string);
    for (int i = 0; i < numBasicFlags ; i++) {
        std::string option = basicFlags[i];
        if (!hasOption(arguments, option)) {
            std::string msg("Please provide the basic arguments;");
            msg += " " + option;
            std::cout << msg << std::endl;
            return 1;
        }
    }
    
    std::string prob_dpath(valueOf(arguments, "-i"));
    std::string appr_name(valueOf(arguments, "-a"));
    std::string appr_dpath = prob_dpath.substr(0, prob_dpath.find_last_of("/"));
    appr_dpath += "/" + appr_name;
    system(("mkdir -p "+ appr_dpath).c_str());
    //
    std::string appr_name_base = appr_name.substr(0, appr_name.find_first_of("-"));
    std::vector<std::string> probFileNames = read_directory(prob_dpath, ".json");
    for (std::string fn: probFileNames) {
        std::string prob_fpath(prob_dpath + "/" + fn);
        Problem *prob = Problem::read_json(prob_fpath);
        std::string postfix = prob->problemName;
        FilePathOrganizer fpo(appr_dpath, postfix);
        std::ifstream is;
        is.open(fpo.solPathCSV);
        if (!is.fail()) {
            is.close();
            continue;
        }
        if (!hasOption(arguments, "-l")) {
            fpo.logPath = "";
        }
        TimeTracker tt;
        std::cout << tt.get_curTime();
        std::cout << "\t" << appr_name << "; " << postfix << std::endl;
        double gap;
        if (appr_name_base == "LP" || appr_name_base == "ILP") {
            bool isTightenModel;
            std::size_t pos = appr_name.find("-");
            if (pos == std::string::npos || appr_name.substr(pos + 1) == "N") {
                isTightenModel = false;
            } else {
                assert(appr_name.substr(pos + 1) == "T");
                isTightenModel = true;
            }
            if (appr_name_base == "LP") {
                LP lp(prob, fpo.logPath, isTightenModel);
                lp.cplex->solve();
                write_solution(prob, fpo, tt, &lp);
            } else {
                ILP ilp(prob, fpo.logPath, isTightenModel);
                ilp.cplex->solve();
                if (ilp.cplex->getStatus() == IloAlgorithm::Infeasible) {
                    ilp.cplex->exportModel(fpo.lpPath.c_str());
                    continue;
                }
                gap = ilp.cplex->getMIPRelativeGap();
                write_solution(prob, fpo, tt, &ilp, gap);
            }
        } else {
            std::vector<std::string> tokens = parseWithDelimiter(appr_name, "-");
            std::vector<std::string> cut_names;
            for (int i = 1; i < tokens.size(); i++) {
                cut_names.push_back(tokens[i]);
            }
            if (tokens[0] == "BnC") {
                BnC bnc(prob, fpo.logPath, cut_names);
                bnc.cplex->solve();
                if (bnc.cplex->getStatus() == IloAlgorithm::Infeasible) {
                    bnc.cplex->exportModel(fpo.lpPath.c_str());
                    continue;
                }
                gap = bnc.cplex->getMIPRelativeGap();
                write_solution(prob, fpo, tt, &bnc, gap);
            } else {
                assert(tokens[0] == "RC");
                RC rc(prob, fpo.logPath, cut_names);
                rc.solve();
            }
        }
    }
    return 0;
}
