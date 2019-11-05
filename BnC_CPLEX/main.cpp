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
#include "GH/IH.hpp"


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

void write_solution(Problem* prob, FilePathOrganizer& fpo, TimeTracker& tt,
                    double objV, double gap, long long int numNodes, int numGenCuts) {
    char row[DEFAULT_BUFFER_SIZE];
    std::fstream fout_csv;
    fout_csv.open(fpo.solPathCSV, std::ios::out);
    fout_csv << "objV,Gap,eliCpuTime,eliWallTime,numNodes,numGenCuts" << "\n";
    sprintf(row, "%f,%f,%f,%f,%lld,%d",
            objV, gap,
            tt.get_elipsedTimeCPU(),
            tt.get_elipsedTimeWall(),
            numNodes, numGenCuts);
    fout_csv << row << "\n";
    fout_csv.close();
}

void write_visitingSeq_arrivalTime(Problem* prob, FilePathOrganizer& fpo,
                                   double** x_ij, double* u_i) {
    std::fstream fout_txt;
    fout_txt.open(fpo.solPathTXT, std::ios::out);
    std::map<int, int> _route;
    for (int i: (*prob).N) {
        for (int j: (*prob).N) {
            if (x_ij[i][j] > 0.5) {
                _route[i] = j;
//                std::cout << "(" << i << "," << j << ")" << std::endl;
            }
        }
    }
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
}

int main(int argc, const char * argv[]) {
    std::vector<std::string> arguments;
    for (int i = 0; i < argc; i++) {
        arguments.push_back(argv[i]);
    }
    std::string basicFlags[] = {"-i", "-a"};
    int numBasicFlags = sizeof(basicFlags) / sizeof(std::string);
    for (int i = 0; i < numBasicFlags; i++) {
        std::string option = basicFlags[i];
        if (!hasOption(arguments, option)) {
            std::string msg("Please provide the basic arguments;");
            msg += " " + option;
            std::cout << msg << std::endl;
            return 1;
        }
    }
    bool enforcementMode = false;
    if (hasOption(arguments, "-ef")) {
        enforcementMode = true;
    }
    std::string prob_dpath(valueOf(arguments, "-i"));
    std::string appr_name(valueOf(arguments, "-a"));
    std::string appr_dpath = prob_dpath.substr(0, prob_dpath.find_last_of("/"));
    appr_dpath += "/" + appr_name;
    system(("mkdir -p " + appr_dpath).c_str());
    //
    std::string appr_name_base = appr_name.substr(0, appr_name.find_first_of("-"));
    std::vector<std::string> probFileNames = read_directory(prob_dpath, ".json");
    for (std::string fn: probFileNames) {
        std::string prob_fpath(prob_dpath + "/" + fn);
        Problem *prob = Problem::read_json(prob_fpath);
        std::string postfix = prob->problemName;
        FilePathOrganizer fpo(appr_dpath, postfix);
        if (!enforcementMode) {
            std::ifstream is;
            std::string handledFiles[] = {fpo.solPathCSV, fpo.lpPath};
            int numAssoFiles = sizeof(handledFiles) / sizeof(std::string);
            bool isHandled = false;
            for (int i = 0; i < numAssoFiles ; i++) {
                is.open(handledFiles[i]);
                if (!is.fail()) {
                    is.close();
                    isHandled = true;
                    break;
                }
            }
            if (isHandled) {
                continue;
            }
        }
        if (!hasOption(arguments, "-l")) {
            fpo.logPath = "";
        }
        TimeTracker tt;
        std::cout << tt.get_curTime();
        std::cout << "\t" << appr_name << "; " << postfix << std::endl;
        if (appr_name_base == "GH") {
            InsertionHeuristic gh(prob, fpo.logPath);
            gh.run();
            double objV = gh.get_objV();
            double gap = -1.0;
            long long int numNodes = 1;
            int numGenCuts = -1;
            write_solution(prob, fpo, tt, objV, gap, numNodes, numGenCuts);
            write_visitingSeq_arrivalTime(prob, fpo, gh.x_ij, gh.u_i);
        } else {
            BaseMM* mm;
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
                    mm = new LP(prob, fpo.logPath, isTightenModel);
                } else {
                    mm = new ILP(prob, fpo.logPath, isTightenModel);
                }
            } else {
                std::vector<std::string> tokens = parseWithDelimiter(appr_name, "-");
                std::vector<std::string> cut_names;
                for (int i = 1; i < tokens.size(); i++) {
                    cut_names.push_back(tokens[i]);
                }
                IloCplex::CutManagement cutManagerType = IloCplex::UseCutForce;
                IloBool isLocalCutAdd = IloFalse;
                if (hasOption(arguments, "-c")) {
                    std::string cutPrmts = valueOf(arguments, "-c");
                    if (cutPrmts.substr(1, 2) == "l") {
                        isLocalCutAdd = IloTrue;
                    } else {
                        assert(cutPrmts.substr(1, 2) == "g");
                    }
                    if (cutPrmts.substr(0, 1) == "i") {
                        cutManagerType = IloCplex::UseCutFilter;
                    } else if(cutPrmts.substr(0, 1) == "u") {
                        cutManagerType = IloCplex::UseCutPurge;
                    } else {
                        assert(cutPrmts.substr(0, 1) == "o");
                    }
                }
                if (tokens[0] == "BnC") {
                    mm = new BnC(prob, fpo.logPath, cut_names, cutManagerType, isLocalCutAdd);
                } else {
                    assert(tokens[0] == "RC");
                    mm = new RC(prob, fpo.logPath, cut_names);
                }
            }
            if (hasOption(arguments, "-t")) {
                int timeLimitSec = std::stoi(valueOf(arguments, "-t"));
                mm->cplex->setParam(IloCplex::TiLim, timeLimitSec);
            }
            if (hasOption(arguments, "-g")) {
                InsertionHeuristic gh(prob, fpo.logPath);
                gh.run();
                if (gh.get_objV() == 0) {
                    double objV = gh.get_objV();
                    double gap = -1.0;
                    long long int numNodes = -1;
                    int numGenCuts = -1;
                    write_solution(prob, fpo, tt, objV, gap, numNodes, numGenCuts);
                    write_visitingSeq_arrivalTime(prob, fpo, gh.x_ij, gh.u_i);
                    continue;
                }
                mm->start_fromGHSol(gh.x_ij, gh.u_i);
            }
            mm->cplex->solve();
            if (mm->cplex->getStatus() == IloAlgorithm::Infeasible) {
                mm->cplex->exportModel(fpo.lpPath.c_str());
                continue;
            }
            if (appr_name_base == "LP" || appr_name_base == "RC") {
                write_solution(prob, fpo, tt, mm);
            } else {
                try {
                    double objV = mm->cplex->getObjValue();
                    double gap = mm->cplex->getMIPRelativeGap();
                    long long int numNodes = mm->cplex->getNnodes64();
                    int numGenCuts = mm->getNumGenCuts();
                    write_solution(prob, fpo, tt, objV, gap, numNodes, numGenCuts);
                    //
                    double **x_ij = new double *[(*prob).N.size()];
                    double *u_i = new double[(*prob).N.size()];
                    for (int i: (*prob).N) {
                        x_ij[i] = new double[(*prob).N.size()];
                    }
                    mm->get_x_ij(x_ij);
                    mm->get_u_i(u_i);
                    write_visitingSeq_arrivalTime(prob, fpo, x_ij, u_i);
                    for (int i: (*prob).N) {
                        delete [] x_ij[i];
                    }
                    delete [] x_ij;
                    delete [] u_i;
                } catch (IloCplex::Exception e) {
                    std::cout << "no incumbent until the time limit" << std::endl;
                    std::fstream fout;
                    fout.open(fpo.lpPath, std::ios::out);
                    fout << "\t No incumbent until the time limit, " << std::stoi(valueOf(arguments, "-t")) << "seconds" << "\n";
                    fout.close();
                }
            }
            delete mm;
        }
//        delete prob;
    }
    return 0;
}
