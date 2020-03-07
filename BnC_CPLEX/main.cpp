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

#include "include/ck_route/Problem.hpp"
#include "include/ck_route/RouteMM.hpp"
#include "src/MathematicalModel/DF.hpp"
#include "src/GH/IH.hpp"
//
#include "ck_util/util.hpp"

#define DEFAULT_BUFFER_SIZE 2048

using namespace rmm;

template<typename Base, typename T>
inline bool instanceof(const T*) {
   return std::is_base_of<Base, T>::value;
}

void write_solution(Problem *prob, FilePathOrganizer &fpo, TimeTracker &tt, RouteMM *mm) {
    char row[DEFAULT_BUFFER_SIZE];
    std::fstream fout_csv;
    fout_csv.open(fpo.solPathCSV, std::ios::out);
    fout_csv << "objV,Gap,eliCpuTime,eliWallTime,notes" << "\n";
    if (instanceof<LP>(mm)) {
        fout_csv << "objV,eliCpuTime,eliWallTime" << "\n";
        sprintf(row, "%f,%f,%f,\"{numRows: %ld,numCols: %ld}\"",
                mm->cplex->getObjValue(),
                tt.get_elapsedTimeCPU(),
                tt.get_elapsedTimeWall(),
                mm->cplex->getNrows(),
                mm->cplex->getNcols());
    } else {
        sprintf(row,
                "%f,%f,%f,%f,\"{\'numNodes\': %lld, \'numGenCuts\': %d, \'time4Sep\': %f,\'time4FDV\': %f, \'num4Sep\': %d}\"",
                mm->cplex->getObjValue(),
                mm->cplex->getMIPRelativeGap(),
                tt.get_elapsedTimeCPU(),
                tt.get_elapsedTimeWall(),
                mm->cplex->getNnodes64(),
                mm->getNumGenCuts(),
                mm->getTime4Sep(),
                mm->getTime4FDV(),
                mm->getNum4Sep());
    }
    fout_csv << row << "\n";
    fout_csv.close();
}

void write_solution(Problem *prob, FilePathOrganizer &fpo, TimeTracker &tt, DF *mm) {
    char row[DEFAULT_BUFFER_SIZE];
    std::fstream fout_csv;
    fout_csv.open(fpo.solPathCSV, std::ios::out);
    fout_csv << "objV,Gap,eliCpuTime,eliWallTime,notes" << "\n";
    sprintf(row,
            "%f,%f,%f,%f",
            mm->cplex->getObjValue(),
            mm->cplex->getMIPRelativeGap(),
            tt.get_elapsedTimeCPU(),
            tt.get_elapsedTimeWall());
    fout_csv << row << "\n";
    fout_csv.close();
}


void write_solution(Problem* prob, FilePathOrganizer& fpo, TimeTracker& tt,
                    double objV, double gap, long long int numNodes, int numGenCuts, double time4Sep) {
    char row[DEFAULT_BUFFER_SIZE];
    std::fstream fout_csv;
    fout_csv.open(fpo.solPathCSV, std::ios::out);
    fout_csv << "objV,Gap,eliCpuTime,eliWallTime,notes" << "\n";
    sprintf(row,
            "%f,%f,%f,%f,\"{numNodes: %lld, numGenCuts: %d, time4Sep: %f}\"",
            objV, gap,
            tt.get_elapsedTimeCPU(),
            tt.get_elapsedTimeWall(),
            numNodes, numGenCuts, time4Sep);
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
        Problem *prob;
        if (appr_name_base != "DF") {
            prob = Problem::read_json(prob_fpath);
        } else {
            prob = Problem::read_json_DF(prob_fpath);
        }
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
            double time4Sep = -1.0;
            write_solution(prob, fpo, tt, objV, gap, numNodes, numGenCuts, time4Sep);
            write_visitingSeq_arrivalTime(prob, fpo, gh.x_ij, gh.u_i);
        } else if (appr_name_base == "DF"){
            DF* mm;
            mm = new DF(prob, fpo.logPath);
            if (hasOption(arguments, "-t")) {
                int timeLimitSec = std::stoi(valueOf(arguments, "-t"));
                mm->cplex->setParam(IloCplex::TiLim, timeLimitSec);
            }
            mm->cplex->solve();
            if (mm->cplex->getStatus() == IloAlgorithm::Infeasible) {
                mm->cplex->exportModel(fpo.lpPath.c_str());
                continue;
            }
            try {
                write_solution(prob, fpo, tt, mm);
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
            delete mm;
        } else {
            RouteMM* mm;
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
                bool isTightenModel = false;
                if (hasOption(arguments, "-ed")) {
                    isTightenModel = true;
                }
                if (tokens[0] == "BnC") {
                    mm = new BnC(prob, fpo.logPath, cut_names, isTightenModel, cutManagerType, isLocalCutAdd, &tt);
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
                    double time4Sep = -1.0;
                    write_solution(prob, fpo, tt, objV, gap, numNodes, numGenCuts, time4Sep);
                    write_visitingSeq_arrivalTime(prob, fpo, gh.x_ij, gh.u_i);
                    continue;
                }
                mm->set_initSol(gh.x_ij, gh.u_i);
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
                    write_solution(prob, fpo, tt, mm);
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
