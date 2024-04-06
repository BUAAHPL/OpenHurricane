/*!
 * \file controllerStringList.cpp
 * \brief Main subroutines for controller string list.
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

#include "controllerStringList.hpp"
#include "controller.hpp"

OpenHurricane::controllerStringList::controllerStringList(const controller &cont) : cont_(cont) {}

hur_nodiscard OpenHurricane::stdStringList
OpenHurricane::controllerStringList::operator()(const string &w) const {
    if (cont_.found(w)) {
        std::string rnl = cont_.findText(w);
        replaceAllMarks(rnl, "\n", " ");
        stdStringList rll;
        if (!rnl.empty()) {
            split(rnl, rll, ",");
            for (integer i = 0; i < rll.size(); ++i) {
                trim(rll[i]);
            }
        }

        return rll;
    }
    return stdStringList();
}
