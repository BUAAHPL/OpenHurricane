/*!
 * \file atomicWeights.cpp
 * \brief Main subroutines for the atomic weigths.
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
#include "atomicWeights.hpp"

namespace OpenHurricane {

    hur_nodiscard integer atomicWeight::checkNum(const char *name) {
        for (integer i = 0; i < nElements_; i++) {
            if (mapName(atomicWeights_[i].name_, name)) {
                return i;
            }
        }
        return -1;
    }

    hur_nodiscard bool atomicWeight::mapName(const char *tabelName, const char *name) {
        integer n = 0;
        while (tabelName[n] != '\0' && name[n] != '\0') {
            if (toupper(tabelName[n]) == toupper(name[n])) {
                n++;
            } else {
                return false;
            }
        }
        return true;
    }

    /*!\brief Search the atomic weight of atomic name[3]
            Return 1 if success
            Return 0 if false.*/
    integer atomicWeight::search(const char *name, real &weight) {
        integer n;
        n = checkNum(name);
        if (n != -1) {
            weight = atomicWeights_[n].weight_;
            return 1;
        }
        return 0;
    }

} // namespace OpenHurricane
