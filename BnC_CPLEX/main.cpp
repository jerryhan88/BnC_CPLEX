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

#include "include/ck_route/Other.hpp"
#include "include/ck_route/Router.hpp"
#include "src/MathematicalModel/DF.hpp"
//
#include "ck_util/util.hpp"

#define DEFAULT_BUFFER_SIZE 2048

using namespace rmm;

template<typename Base, typename T>
inline bool instanceof(const T*) {
   return std::is_base_of<Base, T>::value;
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
        rut::Problem *prob;
        if (appr_name_base != "DF") {
            prob = rut::Problem::read_json(prob_fpath);
        } else {
            prob = rut::Problem::read_json_DF(prob_fpath);
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
        unsigned long time_limit_sec;
        if (hasOption(arguments, "-t")) {
            time_limit_sec = std::stoi(valueOf(arguments, "-t"));
        } else {
            time_limit_sec = ULONG_MAX;
        }
        //
        TimeTracker tt;
        std::cout << tt.get_curTime();
        std::cout << "\t" << appr_name << "; " << postfix << std::endl;
        //
        rut::Solution *sol;
        if (appr_name_base == "GH") {
            rgh::InsertionHeuristic gh(prob, &tt, time_limit_sec, fpo.logPath);
            sol = gh.solve();
        } else if (appr_name_base == "DF"){
            DF* mm;
            mm = new DF(prob, &tt, time_limit_sec, fpo.logPath, fpo.lpPath);
            if (hasOption(arguments, "-t")) {
                int timeLimitSec = std::stoi(valueOf(arguments, "-t"));
                mm->cplex->setParam(IloCplex::TiLim, timeLimitSec);
            }
            sol = mm->solve();
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
                    mm = new LP(prob, &tt, time_limit_sec,
                                fpo.logPath, fpo.lpPath, isTightenModel);
                } else {
                    mm = new ILP(prob, &tt, time_limit_sec,
                                 fpo.logPath, fpo.lpPath, isTightenModel);
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
                    mm = new BnC(prob, &tt, time_limit_sec,
                                 fpo.logPath, fpo.lpPath,
                                 cut_names, isTightenModel,
                                 cutManagerType, isLocalCutAdd);
                } else {
                    assert(tokens[0] == "RC");
                    mm = new RC(prob, &tt, time_limit_sec, fpo.logPath, fpo.lpPath, cut_names);
                }
            }
            if (hasOption(arguments, "-t")) {
                int timeLimitSec = std::stoi(valueOf(arguments, "-t"));
                mm->cplex->setParam(IloCplex::TiLim, timeLimitSec);
            }
            
            if (hasOption(arguments, "-g")) {
                rgh::InsertionHeuristic gh(prob, &tt, time_limit_sec, fpo.logPath);
                sol = gh.solve();
                if (gh.get_objV() == 0) {
                    sol->writeSolCSV(fpo.solPathCSV);
                    sol->writeSolTXT(fpo.solPathTXT);
                    continue;
                }
                mm->set_initSol(gh.x_ij, gh.u_i);
            }
            sol = mm->solve();
            delete mm;
        }
        sol->writeSolCSV(fpo.solPathCSV);
        sol->writeSolTXT(fpo.solPathTXT);
        delete sol;
//        delete prob;
    }
    return 0;
}
