/*!
 * \file clock.cpp
 * \brief Main subroutines for the clock.
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
#include "clock.hpp"
#include "HurFormat.hpp"
#include <iomanip>

const char *OpenHurricane::clock::dayNames[] = {"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"};

hur_nodiscard std::string OpenHurricane::clock::dateTime() {
    time_t timep = getTime();
    struct tm *p = localtime(&timep);
    return hurFormat("%d/%d/%d %s %d:%d:%d", 1900 + p->tm_year, p->tm_mon + 1, p->tm_mday,
                     dayNames[p->tm_wday], p->tm_hour, p->tm_min, p->tm_sec);
}

hur_nodiscard double OpenHurricane::hrClock::elapsedClockTime() const {
    newTime_ = getClockTime();
    std::chrono::duration<double> elapsed = newTime_ - startTime_;
    return elapsed.count();
}

hur_nodiscard double OpenHurricane::hrClock::clockTimeIncrement() const {
    lastTime_ = newTime_;
    newTime_ = getClockTime();
    std::chrono::duration<double> elapsed = newTime_ - lastTime_;
    return elapsed.count();
}