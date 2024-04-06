/*!
 * \file errorAbort.hpp
 * \brief Header of error information output.
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
#pragma once
#include "HurFormat.hpp"
#include "preset.hpp"
namespace OpenHurricane {

    void printError(bool success, char const *const errMsg, char const *const func,
                    const char *const file, const int line);

    void printWarning(char const *const errMsg, char const *const func, const char *const file,
                      const int line);

    void printLogInfo(char const *const errMsg, char const *const func, const char *const file,
                      const int line);

    void printPLogInfo(char const *const errMsg, char const *const func, const char *const file,
                       const int line);
    void printPInfo(char const *const errMsg, char const *const func, const char *const file,
                    const int line);

    void printLogWarning(char const *const errMsg, char const *const func, const char *const file,
                         const int line);

    void printPLogWarning(char const *const errMsg, char const *const func, const char *const file,
                          const int line);

    void printPWarning(char const *const errMsg, char const *const func, const char *const file,
                       const int line);

    void printLogDebug(char const *const errMsg, char const *const func, const char *const file,
                       const int line);
    void printPLogDebug(char const *const errMsg, char const *const func, const char *const file,
                        const int line);
    void printPDebug(char const *const errMsg, char const *const func, const char *const file,
                     const int line);

    void printLogError(char const *const errMsg, char const *const func, const char *const file,
                       const int line);

    void printPLogError(char const *const errMsg, char const *const func, const char *const file,
                        const int line);

    void printPError(char const *const errMsg, char const *const func, const char *const file,
                     const int line);

    void printfatalError(char const *const errMsg, char const *const func, const char *const file,
                         const int line);

#ifndef checkWarning
#define checkWarning(warningMsg) printWarning((warningMsg), HUR_FUNC, HUR_FILE, HUR_LINE)
#endif // !checkWarning

#ifndef checkWarningStr
#define checkWarningStr(warningMsgStr) \
    printWarning((warningMsgStr.c_str()), HUR_FUNC, HUR_FILE, HUR_LINE)
#endif // !checkWarningStr

#ifndef errorAbort
#define errorAbort(errMsg) printError(false, (errMsg), HUR_FUNC, HUR_FILE, HUR_LINE)
#endif // !errorAbort
#ifndef errorAbortStr
#define errorAbortStr(errMsgStr) \
    printError(false, (errMsgStr.c_str()), HUR_FUNC, HUR_FILE, HUR_LINE)
#endif // !errorAbortStr

#ifdef HUR_NONE_LOGGER
    /**
     * \brief Output information to log file (optional).
     */
#define LInfo(fmt, ...)
    /**
     * \brief Output information to screen.
     */
#define PInfo(fmt, ...) Pout(fmt, ##__VA_ARGS__)
    /**
     * \brief Output information to screen and log file (optional).
     */
#define PLInfo(fmt, ...) Pout(fmt, ##__VA_ARGS__)

    /**
     * \brief Output warning information to log file (optional).
     */
#define LWarning(fmt, ...)

    /**
     * \brief Output warning information to screen.
     */
#define PWarning(fmt, ...) \
    printPWarning(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)
    /**
     * \brief Output warning information to screen and log file (optional).
     */
#define PLWarning(fmt, ...) \
    printPWarning(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)

#define LDebug(fmt, ...)
#define PDebug(fmt, ...)
#define PLDebug(fmt, ...)

#elif defined(HUR_LESS_LOGGER) || defined(HUR_FULL_LOGGER)
    /**
     * \brief Output information to log file (optional).
     */
#define LInfo(fmt, ...) \
    printLogInfo(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)
    /**
     * \brief Output information to screen.
     */
#define PInfo(fmt, ...) Pout(fmt, ##__VA_ARGS__)
    /**
     * \brief Output information to screen and log file (optional).
     */
#define PLInfo(fmt, ...) \
    printPLogInfo(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)

    /**
     * \brief Output warning information to log file (optional).
     */
#define LWarning(fmt, ...) \
    printLogWarning(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)
    /**
     * \brief Output warning information to screen.
     */
#define PWarning(fmt, ...) \
    printPWarning(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)
    /**
     * \brief Output warning information to screen and log file (optional).
     */
#define PLWarning(fmt, ...) \
    printPLogWarning(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)
    /**
     * \brief Output debug information to log file (optional).
     */
#define LDebug(fmt, ...) \
    printLogDebug(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)

    /**
     * \brief Output debug information to screen (optional).
     */
#define PDebug(fmt, ...) \
    printPDebug(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)

    /**
     * \brief Output debug information to screen and log file (optional).
     */
#define PLDebug(fmt, ...) \
    printPLogDebug(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)
#endif // HUR_NONE_LOGGER

#if defined(HUR_LESS_LOGGER) || defined(HUR_FULL_LOGGER)
    /**
     * \brief Output error information to log file.
     * And will not abort the program.
     */
#define LError(fmt, ...) \
    printLogError(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)

    /**
     * \brief Output error information to screen.
     * And will not abort the program.
     */
#define PError(fmt, ...) \
    printPError(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)

    /**
     * \brief Output error information to screen and log file.
     * And will not abort the program.
     */
#define PLError(fmt, ...) \
    printPLogError(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)
#else
    /**
     * \brief Output error information to log file.
     * And will not abort the program.
     */
#define LError(fmt, ...)
#define PError(fmt, ...) \
    printPError(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)
#define PLError(fmt, ...) \
    printPError(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)

#endif // HUR_NONE_LOGGER

    /**
     * \brief Output fatal information to screen and log file(optional).
     * And abort the program.
     */
#define LFatal(fmt, ...) \
    printfatalError(hurFormat(fmt, ##__VA_ARGS__).c_str(), HUR_FUNC, HUR_FILE, HUR_LINE)

} // namespace OpenHurricane
