#include "logFile.hpp"
/*!
 * \file logFile.cpp
 * \brief The subroutines and functions of log files
 * \author Rao Sihang
 * \version V2.0.0
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

#include "HurFormat.hpp"
#include "HurMPIBase.hpp"
#include "clock.hpp"
#include "logFile.hpp"
#include "version.hpp"
#include <iomanip>

const OpenHurricane::fileName OpenHurricane::loggerType::defaultFileName("OpenHurricane.log");

bool OpenHurricane::report = false;

OpenHurricane::loggerType OpenHurricane::PLogs;

OpenHurricane::loggerType::loggerType()
    : fos_(IOsstream::streamFormat::ASCII_FORMAT, IOsstream::openOptions::CREATE_MASTER,
           std::ios_base::out) {
    if (!report) {
        fos_.close();
    }
}

void OpenHurricane::loggerType::initializing(fileName fn) {
    if (report) {
        if (!fn.isAbsolute()) {
            fn.toAbsolute();
        }
        fos_.changeFileName(fn);
        fos_.open();
    }
}

void OpenHurricane::loggerType::log(const logTypes::logType lt, char const *const errMsg,
                                char const *const func, const char *const file, const int line) {
    if (report && fos_.isOpen()) {
        fileName curfileN(file);
        // [time][Process ID][collective type][logType][file name][function name : line number]: errmsg
        const char *fmt = "[%s][PC:%d][collective:%s][%s][%s][%s:%4d]: %s";
        std::string str;
        if (HurMPIBase::isCollectiveCall()) {
            if (HurMPIBase::master()) {
                str = hurFormat(fmt, clock::dateTime().c_str(), HurMPIBase::getProcRank(), "Y",
                                type(lt).c_str(), curfileN.name().c_str(), func, line, errMsg);
                fos_.os() << str << std::endl << std::flush;
            }
        } else {
            str = hurFormat(fmt, clock::dateTime().c_str(), HurMPIBase::getProcRank(), "N",
                            type(lt).c_str(), curfileN.name().c_str(), func, line, errMsg);
            fos_.os() << str << std::endl << std::flush;
        }
    }
}

std::string OpenHurricane::loggerType::type(const logTypes::logType lt) {
    switch (lt) {
    case (logTypes::INFOTYPE):
        return "INFO";
        break;
    case (logTypes::DEBUGTYPE):
        return "DEBUG";
        break;
    case (logTypes::WARNINGTYPE):
        return "WARNING";
        break;
    case (logTypes::ERRORTYPE):
        return "ERROR";
        break;
    case (logTypes::FATALTYPE):
        return "FATAL";
        break;
    }
    return " UNKNOWN: ";
}