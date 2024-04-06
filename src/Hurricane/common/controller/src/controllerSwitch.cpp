/*!
 * \file controllerSwitch.cpp
 * \brief Main subroutines for controller switch.
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

#include "controllerSwitch.hpp"
#include "controller.hpp"

OpenHurricane::controllerSwitch::controllerSwitch(const controller &cont) : cont_(cont) {}

hur_nodiscard bool OpenHurricane::controllerSwitch::operator()(const string &w, const bool defaultOption,
                                             const char *on, const char *off) const {
    if (cont_.found(w)) {
        std::string onstr;
        std::string offstr;
        onstr = on;
        offstr = off;
        stringToUpperCase(onstr);
        stringToUpperCase(offstr);

        const auto cpw = cont_.findWord(w);
        string pw = cpw;
        trim(pw);
        stringToUpperCase(pw);

        if (pw == onstr) {
            return true;
        } else if (pw == offstr) {
            return false;
        } else {
            errorAbortStr(("Unknown item: " + cpw + " in controller: " + cont_.name()));
        }
    }
    return defaultOption;
}
