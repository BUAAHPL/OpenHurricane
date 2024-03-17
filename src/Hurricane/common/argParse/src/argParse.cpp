/*!
 * \file argParse.cpp
 * \brief The subroutines and functions of arg parsing
 * \author Rao Sihang
 * \version V2.0.0
 * \date 2022.05.02
 *
 * OpenHurricane: Open parts of Hurricane project (Highly Universal Rocket & Ramjet sImulation Codes for ANalysis and Evaluation)
 * \copyright Copyright (C) 2019-2023, Prof. Xu Xu's group at Beihang University.
 *
 * License
 *		This file is part of OpenHurricane
 *
 *		OpenHurricane is free software: you can redistribute it and/or modify it
 *		under the terms of the GNU General Public License as published by
 *		the Free Software Foundation, either version 3 of the License, or
 *		(at your option) any later version.
 *
 *		OpenHurricane is distributed in the hope that it will be useful, but WITHOUT
 *		ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *		FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *		for more details.
 *
 *		You should have received a copy of the GNU General Public License
 *		along with OpenHurricane.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */
#include "argParse.hpp"
#include "Pout.hpp"
#include "logFile.hpp"
#include "reordering.hpp"
#include "spaceSchemeOptions.hpp"
#include "version.hpp"
#include <cctype>

OpenHurricane::string OpenHurricane::argParse::executable_;
OpenHurricane::fileName OpenHurricane::argParse::caseFileName_;
OpenHurricane::fileName OpenHurricane::argParse::dataFileName_;
OpenHurricane::fileName OpenHurricane::argParse::outputFileName_;
OpenHurricane::fileName OpenHurricane::argParse::residualFileName_;
OpenHurricane::fileName OpenHurricane::argParse::executePath_;
OpenHurricane::fileName OpenHurricane::argParse::inputPath_;
OpenHurricane::fileName OpenHurricane::argParse::outputPath_;
bool OpenHurricane::argParse::isASCII_(false);
bool OpenHurricane::argParse::isGUIMode_(false);
bool OpenHurricane::argParse::hasGPU_(false);
int OpenHurricane::argParse::nGPUPerMachine_(0);

std::map<std::string, OpenHurricane::argOptions::options> OpenHurricane::argParse::validOptions;

bool OpenHurricane::argParse::gpuZeroCopy_(false);
bool OpenHurricane::argParse::tecplotFileWriteOnMaster_(false);
bool OpenHurricane::argParse::noWriteAtEnd_(false);
bool OpenHurricane::argParse::isSkewCorrect_(false);
bool OpenHurricane::argParse::checkWriteSolution_(false);

OpenHurricane::argParse::initValidTables::initValidTables() {
    argParse::addValidOptions("HELP", OpenHurricane::argOptions::HELP);
    argParse::addValidOptions("H", OpenHurricane::argOptions::HELP);
    argParse::addValidOptions("DATAFORMAT", OpenHurricane::argOptions::DATAFORMAT);
    argParse::addValidOptions("DF", OpenHurricane::argOptions::DATAFORMAT);
    argParse::addValidOptions("L", OpenHurricane::argOptions::LOGFILE);
    argParse::addValidOptions("LOG", OpenHurricane::argOptions::LOGFILE);
    argParse::addValidOptions("LOGFILE", OpenHurricane::argOptions::LOGFILE);
    argParse::addValidOptions("V", OpenHurricane::argOptions::VERSION);
    argParse::addValidOptions("VERSION", OpenHurricane::argOptions::VERSION);
    argParse::addValidOptions("C", OpenHurricane::argOptions::CASE);
    argParse::addValidOptions("CASE", OpenHurricane::argOptions::CASE);
    argParse::addValidOptions("D", OpenHurricane::argOptions::DATA);
    argParse::addValidOptions("DATA", OpenHurricane::argOptions::DATA);
    argParse::addValidOptions("CD", OpenHurricane::argOptions::CD);
    argParse::addValidOptions("I", OpenHurricane::argOptions::INPUT);
    argParse::addValidOptions("INPUT", OpenHurricane::argOptions::INPUT);
    argParse::addValidOptions("IP", OpenHurricane::argOptions::INPUTPATH);
    argParse::addValidOptions("O", OpenHurricane::argOptions::OUT);
    argParse::addValidOptions("OUT", OpenHurricane::argOptions::OUT);
    argParse::addValidOptions("INPUTPATH", OpenHurricane::argOptions::INPUTPATH);
    argParse::addValidOptions("OUTPUTPATH", OpenHurricane::argOptions::OUTPUTPATH);
    argParse::addValidOptions("OP", OpenHurricane::argOptions::OUTPUTPATH);
    argParse::addValidOptions("R", OpenHurricane::argOptions::RESIDUAL);
    argParse::addValidOptions("RESIDUAL", OpenHurricane::argOptions::RESIDUAL);
    argParse::addValidOptions("GUI", OpenHurricane::argOptions::GUI);
    argParse::addValidOptions("NGPU", OpenHurricane::argOptions::nGPU);
    argParse::addValidOptions("GPUZEROCOPY", OpenHurricane::argOptions::GPUZEROCOPY);
    argParse::addValidOptions("TECPLOTFILEWRITEONMASTER",
                              OpenHurricane::argOptions::TECPLOTFILEWRITEONMASTER);
    argParse::addValidOptions("NOWRITEATEND", OpenHurricane::argOptions::noWRITEAtEND);
    argParse::addValidOptions("SKEWNESSCORRECTION", OpenHurricane::argOptions::SKEWNESSCORRECTION);
    argParse::addValidOptions("CHECKWRITESOLUTION", OpenHurricane::argOptions::CHECKWRITESOLUTION);
}

OpenHurricane::argParse::initValidTables dummyInitValidTables;

void OpenHurricane::argParse::addValidOptions(const std::string &opt, const argOptions::options param) {
    validOptions.emplace(opt, param);
}

OpenHurricane::argParse::argParse(int &argc, char **&argv) : args_(argc) {
    defaultName();
    for (int argI = 0; argI < argc; argI++) {
        args_[argI] = argv[argI];
    }
    args_[0] = fileName(argv[0]);
    executable_ = fileName(argv[0]).name();

    executePath_ = getCurrentPath();
    inputPath_ = executePath_;
    outputPath_ = executePath_;
    argListStr_ = args_[0];

    for (int argI = 1; argI < args_.size(); argI++) {
        argListStr_ += ' ';
        argListStr_ += args_[argI];
        if (args_[argI][0] == '-') {
            const char *optionName = &args_[argI][1];
            std::string optionNameStr(optionName);
            stringToUpperCase(optionNameStr);
            std::map<std::string, argOptions::options>::iterator iter;
            iter = validOptions.find(optionNameStr);
            if ((iter != validOptions.end()) && (!optionNameStr.empty())) {
                if (argI + 1 >= args_.size()) {
                    options_.emplace(iter->second, "");
                    break;
                }

                if (args_[argI + 1][0] == '-') {
                    options_.emplace(iter->second, "");
                    continue;
                }

                ++argI;
                argListStr_ += ' ';
                argListStr_ += args_[argI];
                options_.emplace(iter->second, args_[argI]);
            } else {
                PWarning("    Warning: unknown argument: -%s\n", optionName);
                options_.emplace(argOptions::NON, "");
            }
        } else {
            PWarning("    Warning: It's not an argument: %s\n", args_[argI].c_str());
        }
    }

    parse();
}

OpenHurricane::argParse::argParse(int &argc, char **&argv, const char *prefix) : args_(argc) {
    defaultName(prefix);

    for (int argI = 0; argI < argc; argI++) {
        args_[argI] = argv[argI];
    }
    args_[0] = fileName(argv[0]);
    executable_ = fileName(argv[0]).name();
    executePath_ = getCurrentPath();
    inputPath_ = executePath_;
    outputPath_ = executePath_;
    argListStr_ = args_[0];

    for (int argI = 1; argI < args_.size(); argI++) {
        argListStr_ += ' ';
        argListStr_ += args_[argI];
        if (args_[argI][0] == '-') {
            const char *optionName = &args_[argI][1];
            std::string optionNameStr(optionName);
            stringToUpperCase(optionNameStr);
            std::map<std::string, argOptions::options>::iterator iter;
            iter = validOptions.find(optionNameStr);
            if ((iter != validOptions.end()) && (!optionNameStr.empty())) {
                if (argI + 1 >= args_.size()) {
                    options_.emplace(iter->second, "");
                    break;
                }

                if (args_[argI + 1][0] == '-') {
                    options_.emplace(iter->second, "");
                    continue;
                }

                ++argI;
                argListStr_ += ' ';
                argListStr_ += args_[argI];
                options_.emplace(iter->second, args_[argI]);
            } else {
                PWarning("  Warning: unknown argument: -%s\n", optionName);
                options_.emplace(argOptions::NON, "");
            }
        } else {
            PWarning("  Warning: It's not an argument: %s\n", args_[argI].c_str());
        }
    }

    parse();
}

OpenHurricane::argParse::argParse(int &argc, char **&argv, const string &prefix) : args_(argc) {
    defaultName(prefix);

    for (int argI = 0; argI < argc; argI++) {
        args_[argI] = argv[argI];
    }
    args_[0] = fileName(argv[0]);
    executable_ = fileName(argv[0]).name();
    executePath_ = getCurrentPath();
    inputPath_ = executePath_;
    outputPath_ = executePath_;
    argListStr_ = args_[0];

    for (int argI = 1; argI < args_.size(); argI++) {
        argListStr_ += ' ';
        argListStr_ += args_[argI];
        if (args_[argI][0] == '-') {
            const char *optionName = &args_[argI][1];
            std::string optionNameStr(optionName);
            stringToUpperCase(optionNameStr);
            std::map<std::string, argOptions::options>::iterator iter;
            iter = validOptions.find(optionNameStr);
            if ((iter != validOptions.end()) && (!optionNameStr.empty())) {
                if (argI + 1 >= args_.size()) {
                    options_.emplace(iter->second, "");
                    break;
                }

                if (args_[argI + 1][0] == '-') {
                    options_.emplace(iter->second, "");
                    continue;
                }

                ++argI;
                argListStr_ += ' ';
                argListStr_ += args_[argI];
                options_.emplace(iter->second, args_[argI]);
            } else {
                PWarning("  Warning: unknown argument: -%s\n", optionName);
                options_.emplace(argOptions::NON, "");
            }
        } else {
            PWarning("  Warning: It's not an argument: %s\n", args_[argI].c_str());
        }
    }

    parse();
}

OpenHurricane::argParse::~argParse() noexcept {
    clear();
}

void OpenHurricane::argParse::parse() {
    mapType::iterator iter;
    iter = options_.find(argOptions::HELP);
    if (iter != options_.end()) {
        printHelp();
        HurMPI::finalize();
        exit(EXIT_SUCCESS);
    }

    iter = options_.find(argOptions::VERSION);
    if (iter != options_.end()) {
        printVersion();
        HurMPI::finalize();
        exit(EXIT_SUCCESS);
    }

    iter = options_.find(argOptions::DATAFORMAT);
    if (iter != options_.end()) {
        const auto df = iter->second;
        if (df == "0") {
            isASCII_ = true;
        } else if (df != "1") {
            PWarning("Invalid parameter for data format setting. To set the data format output by "
                     "this executable program.\n"
                     "              0 - ASCII;\n"
                     "              1 - Binary.\n"
                     "          If not set, the data would be output in binary.");
        }
    }

    iter = options_.find(argOptions::INPUTPATH);
    if (iter != options_.end()) {
        fileName fn = iter->second;
        if (!fn.empty()) {
            inputPath_ = fn;
            outputPath_ = inputPath_;
        } else {
            LFatal("Please specify the input path");
        }
    }

    iter = options_.find(argOptions::INPUT);
    if (iter != options_.end()) {
        fileName fn = iter->second;
        if (!fn.empty()) {
            if (fn.isAbsolute()) {
                caseFileName_ = fn;
                inputPath_ = caseFileName_.parentPath();
                outputPath_ = inputPath_;
            } else {
                caseFileName_ = inputPath_ / fn;
                inputPath_ = caseFileName_.parentPath();
                outputPath_ = inputPath_;
            }
            dataFileName_ = fn.name(true);
            outputFileName_ = dataFileName_;
            dataFileName_ = dataFileName_ + ".h5";
            outputFileName_ = outputFileName_ + "out.dat";
        } else {
            LFatal("Please specify the case file");
        }
    } else {
        caseFileName_ = inputPath_ / caseFileName_;
    }

    iter = options_.find(argOptions::OUTPUTPATH);
    if (iter != options_.end()) {
        fileName fn = iter->second;
        if (!fn.empty()) {
            outputPath_ = fn;
        } else {
            LFatal("Please specify the input path");
        }
    }

    iter = options_.find(argOptions::OUT);
    if (iter != options_.end()) {
        fileName fn = iter->second;
        if (!fn.empty()) {
            if (fn.isAbsolute()) {
                dataFileName_ = fn;
                outputPath_ = dataFileName_.parentPath();
            } else {
                dataFileName_ = outputPath_ / fn;
            }
            outputFileName_ = fn.name(true);
            outputFileName_ += "out.dat";
            outputFileName_ = outputPath_ / outputFileName_;
        } else {
            LFatal("Please specify the case file");
        }
    } else {
        dataFileName_ = outputPath_ / dataFileName_;
        outputFileName_ = outputPath_ / outputFileName_;
    }

    iter = options_.find(argOptions::RESIDUAL);
    if (iter != options_.end()) {
        fileName fn = iter->second;
        if (!fn.empty()) {
            if (fn.isAbsolute()) {
                residualFileName_ = fn;
            } else {
                residualFileName_ = outputPath_ / fn;
            }
        } else {
            LFatal("Please specify the residual file");
        }
    }

    iter = options_.find(argOptions::LOGFILE);
    if (iter != options_.end()) {
        report = true;
        fileName fn = iter->second;
        if (!fn.empty()) {
            PLogs.initializing(fn);
        } else {
            PLogs.initializing();
        }
    }

    iter = options_.find(argOptions::GUI);
    if (iter != options_.end()) {
        isGUIMode_ = true;
    }

    iter = options_.find(argOptions::nGPU);
    if (iter != options_.end()) {
        hasGPU_ = true;
        auto num = iter->second;
        if (!num.empty()) {
            std::stringstream sstr(num);
            sstr >> nGPUPerMachine_;
        } else {
            nGPUPerMachine_ = 0;
        }
    }

    iter = options_.find(argOptions::GPUZEROCOPY);
    if (iter != options_.end()) {
        gpuZeroCopy_ = true;
    }

    iter = options_.find(argOptions::TECPLOTFILEWRITEONMASTER);
    if (iter != options_.end()) {
        tecplotFileWriteOnMaster_ = true;
    }

    iter = options_.find(argOptions::noWRITEAtEND);
    if (iter != options_.end()) {
        noWriteAtEnd_ = true;
    }

    iter = options_.find(argOptions::SKEWNESSCORRECTION);
    if (iter != options_.end()) {
        isSkewCorrect_ = true;
    }

    iter = options_.find(argOptions::CHECKWRITESOLUTION);
    if (iter != options_.end()) {
        checkWriteSolution_ = true;
    }
}

void OpenHurricane::argParse::printHelp() const {
    Pout("\nUsage: ./%s [opts] <opts args> ...\n\n", executable_.c_str());

    Pout("Options\n");
    Pout("========================================\n");

    auto printHp = [](const char *args, const char *msg, const char *abbre = nullptr,
                      const char *option = nullptr, const int fixWidth = 82) {
        const auto len = static_cast<int>(strlen(msg));
        const std::string fixWidthStr = "          ";
        std::string msgStr;
        auto ls = len / fixWidth;
        if (ls <= 0) {
            msgStr += fixWidthStr;
            msgStr += msg;
        } else {
            int cStart = 0;
            for (int i = 0; i <= ls; ++i) {
                auto cSize = fixWidth;
                if (i == ls) {
                    cSize = len - cStart;
                }

                char *c = new char[cSize + 1]();
                strncpy(c, msg + cStart, cSize);
                c[cSize] = '\0';
                msgStr += fixWidthStr;
                msgStr += c;
                if (i != ls) {
                    if (*(msg + cStart + cSize) == '.') {
                        msgStr += ".";
                        cSize++;
                    } else if (*(msg + cStart + cSize - 1) != ' ' &&
                               *(msg + cStart + cSize) != ' ') {
                        msgStr += "-";
                    }
                    msgStr += "\n";
                }
                HurDeleteDynArray(c);
                cStart += cSize;
            }
        }

        if (abbre == nullptr) {
            Pout(" -%s\n%s\n\n", args, msgStr.c_str());
        } else if (option == nullptr) {
            Pout(" -%s    (-%s)\n%s\n\n", args, abbre, msgStr.c_str());
        } else {
            Pout(" -%s    (-%s)   <%s>\n%s\n\n", args, abbre, option, msgStr.c_str());
        }
    };
    printHp("help", "Print this help information on this executable program.", "h");
    printHp("input",
            hurFormat("Input the case file specified by <file name>. If not specify the "
                      "case file name, the default case filename: \"%s\" would be used.",
                      caseFileName_.c_str())
                .c_str(),
            "I", "file name");

    printHp("inputpath",
            "Set the file path specified by <file path>. If not "
            "specify the file path,the case file and executable file should be included in "
            "the same directory.",
            "IP", "file path");
    printHp(
        "output",
        hurFormat("Output the data file specified by <file name>. If not "
                  "specify the case file name, the default case filename: \"%s\" would be used.",
                  dataFileName_.c_str())
            .c_str(),
        "O", "file name");
    printHp("outputpath",
            "Set the file path specified by <file path>. If not "
            "specify the file path,the data file and executable file should be included in "
            "the same directory.",
            "OP", "file path");
    printHp("log",
            hurFormat("Turn on the log file system on this executable program, and use the "
                      "default log file: \"%s\" if the file name is not specified.",
                      loggerType::defaultFileName.c_str())
                .c_str(),
            "L", "file name (optional)");
    printHp("dataformat",
            "To set the data format output by this executable program: 0 - ASCII; 1 - Binary. If "
            "not set, the data would be output in binary.",
            "DF", "0 | 1");
    printHp("version", "Print version information on this executable program.", "V");
    printHp("nGPU", "To set the number of GPUs per machine.");
    printHp("gpuZeroCopy",
            "To set on zero copy from cpu to gpu. (Map cpu hostmemory to CUDA devices)");
    printHp("tecplotFileWriteOnMaster",
            "Set the tecplot file writting on master process. Default is false.");
    printHp("noWriteAtEnd", "To specify do not write the output and relay files at the "
                                        "end of the computation. Default is false.");
    printHp("skewnessCorrection",
            "To specify the implementation of a correction for the interpolated value at the face "
            "centre when the face is skewed. Default is false.");

    printHp("checkWriteSolution",
            "To specify checking the solution for writting out is NaN or INF?. Default is false.");
    printVersion();
}

void OpenHurricane::argParse::printVersion() const {
    if (HurMPI::master()) {
        programName::printVersion();
    }
    HurMPI::barrierDefaultComm();
}

void OpenHurricane::argParse::clear() noexcept {
    args_.clear();
    options_.clear();
    argListStr_.clear();
}