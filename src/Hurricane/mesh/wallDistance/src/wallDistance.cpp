/*!
 * \file wallDistance.cpp
 * \brief Main subroutines for search procedures of wall distance.
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
#include "wallDistance.hpp"

OpenHurricane::wallDistance::wallDistance(const runtimeMesh &mesh, const controller &cont,
                                      faceBCType::bcTypes tp)
    : mesh_(mesh), zoneIdList_(), zoneType_(tp),
      wallDist_(object("wallDist", mesh, object::WRITE_OUTPUT), mesh), dmPtr_(nullptr) {
    const auto &fzl = mesh.faceZones();
    integer countW = Zero;
    for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
        if (fzl[fzi].bcType() == tp) {
            zoneIdList_.append(fzi);
            countW++;
        }
    }

    if (countW == 0) {
        Pout << std::endl
             << "    Info: there is no wall face zone, and the neareast wall "
                "distance is set to a great value. "
             << std::endl
             << std::endl;
        wallDist_ = large;
        return;
    }

    Pout << std::endl << "    Info: getting distance from cells to nearest wall. " << std::endl;
    dmPtr_ = distanceMethod::creator(cont, mesh, zoneIdList_);

    dmPtr_->getDistance(wallDist_);
}
