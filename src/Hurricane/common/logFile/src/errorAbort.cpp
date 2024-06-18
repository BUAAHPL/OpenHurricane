#include "errorAbort.hpp"
/*!
 * \file errorAbort.cpp
 * \brief Main subroutines for error information output.
 * \author Rao Sihang
 * \version V1.0.0
 * \date 2022.05.02
 *
 * OpenHurricane: Open parts of Hurricane project (Highly Universal Rocket & Ramjet sImulation Codes for ANalysis and Evaluation)
 * \copyright Copyright (C) 2019-2024, Prof. Xu Xu's group at Beihang University.
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

#include "HurMPIBase.hpp"
#include "errorAbort.hpp"
#include "logFile.hpp"

void OpenHurricane::printError(bool success, char const *const errMsg, char const *const func,
                           const char *const file, const int line) {
    if (!success) {
        std::string functionName = func;
#ifdef HUR_DEBUG
        functionName += " of ";
        functionName += file;
        functionName += " at line: ";
        functionName += std::to_string(line);
#endif // HUR_DEBUG
        HurMPIBase::error(errMsg, functionName);
    }
}

void OpenHurricane::printWarning(char const *const errMsg, char const *const func,
                             const char *const file, const int line) {
    std::string functionName = func;
#ifdef HUR_DEBUG
    functionName += " of ";
    functionName += file;
    functionName += " at line: ";
    functionName += std::to_string(line);
#endif // HUR_DEBUG
    HurMPIBase::warning(errMsg, functionName);
}

void OpenHurricane::printLogInfo(char const *const errMsg, char const *const func,
                             const char *const file, const int line) {
    PLogs(logTypes::INFOTYPE, errMsg, func, file, line);
}
void OpenHurricane::printPLogInfo(char const *const errMsg, char const *const func,
                              const char *const file, const int line) {
    PLogs(logTypes::INFOTYPE, errMsg, func, file, line);
    Pout("%s\n", errMsg);
}

void OpenHurricane::printPInfo(char const *const errMsg, char const *const func, const char *const file,
                           const int line) {
    Pout("%s\n", errMsg);
}

void OpenHurricane::printLogWarning(char const *const errMsg, char const *const func,
                                const char *const file, const int line) {
    PLogs(logTypes::WARNINGTYPE, errMsg, func, file, line);
}

void OpenHurricane::printPLogWarning(char const *const errMsg, char const *const func,
                                 const char *const file, const int line) {
    PLogs(logTypes::WARNINGTYPE, errMsg, func, file, line);

#ifdef HUR_FULL_LOGGER
    std::string functionName = func;
#ifdef HUR_DEBUG
    functionName += " of ";
    functionName += file;
    functionName += " at line: ";
    functionName += std::to_string(line);
#endif // HUR_DEBUG
    HurMPIBase::warning(errMsg, functionName);
#elif defined(HUR_LESS_LOGGER)
    Pout("%s in function: %s\n", errMsg, func);
#else
    Pout("%s\n", errMsg);
#endif // HUR_FULL_LOGGER
}

void OpenHurricane::printPWarning(char const *const errMsg, char const *const func,
                              const char *const file, const int line) {
#ifdef HUR_FULL_LOGGER
    std::string functionName = func;
#ifdef HUR_DEBUG
    functionName += " of ";
    functionName += file;
    functionName += " at line: ";
    functionName += std::to_string(line);
#endif // HUR_DEBUG
    HurMPIBase::warning(errMsg, functionName);
#elif defined(HUR_LESS_LOGGER)
    Pout("%s in function: %s\n", errMsg, func);
#else
    Pout("%s\n", errMsg);
#endif // HUR_FULL_LOGGER
}

void OpenHurricane::printLogDebug(char const *const errMsg, char const *const func,
                              const char *const file, const int line) {
    PLogs(logTypes::DEBUGTYPE, errMsg, func, file, line);
}

void OpenHurricane::printPLogDebug(char const *const errMsg, char const *const func,
                               const char *const file, const int line) {
    PLogs(logTypes::DEBUGTYPE, errMsg, func, file, line);
#if defined(HUR_FULL_LOGGER) || defined(HUR_LESS_LOGGER)
    fileName curfileN(file);
    Pout("%s in function %s of file %s in line %d\n", errMsg, func, curfileN.name().c_str(), line);
#else
    Pout("%s\n", errMsg);
#endif // HUR_FULL_LOGGER
}

void OpenHurricane::printPDebug(char const *const errMsg, char const *const func,
                            const char *const file, const int line) {
#if defined(HUR_FULL_LOGGER) || defined(HUR_LESS_LOGGER)
    fileName curfileN(file);
    Pout("%s in function %s of file %s in line %d\n", errMsg, func, curfileN.name().c_str(), line);
#else
    Pout("%s\n", errMsg);
#endif // HUR_FULL_LOGGER
}

void OpenHurricane::printLogError(char const *const errMsg, char const *const func,
                              const char *const file, const int line) {
    PLogs(logTypes::ERRORTYPE, errMsg, func, file, line);
}

void OpenHurricane::printPLogError(char const *const errMsg, char const *const func,
                               const char *const file, const int line) {
    PLogs(logTypes::ERRORTYPE, errMsg, func, file, line);
    std::string functionName = func;
#ifdef HUR_DEBUG
    functionName += " of ";
    functionName += file;
    functionName += " at line: ";
    functionName += std::to_string(line);
#endif // HUR_DEBUG
    HurMPIBase::error(errMsg, functionName, false);
}
void OpenHurricane::printPError(char const *const errMsg, char const *const func,
                            const char *const file, const int line) {
    std::string functionName = func;
#ifdef HUR_DEBUG
    functionName += " of ";
    functionName += file;
    functionName += " at line: ";
    functionName += std::to_string(line);
#endif // HUR_DEBUG
    HurMPIBase::error(errMsg, functionName, false);
}

void OpenHurricane::printfatalError(char const *const errMsg, char const *const func,
                                const char *const file, const int line) {
    if (report) {
        PLogs(logTypes::FATALTYPE, errMsg, func, file, line);
    }
    std::string functionName = func;
#ifdef HUR_DEBUG
    functionName += " of ";
    functionName += file;
    functionName += " at line: ";
    functionName += std::to_string(line);
#endif // HUR_DEBUG
    HurMPIBase::error(errMsg, functionName);
}
