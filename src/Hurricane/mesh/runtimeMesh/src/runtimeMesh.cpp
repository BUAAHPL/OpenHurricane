/*!
 * \file runtimeMesh.cpp
 * \brief Main subroutines for runtime mesh.
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

#include "runtimeMesh.hpp"
#include "NullRefObj.hpp"
#include "boundaries.hpp"
#include "cellMesh.hpp"
#include "faceMesh.hpp"
#include "geometryArrays.hpp"

void OpenHurricane::runtimeMesh::createRegions(const controller &cont) const {
    if (!cont.found("meshSet")) {
        checkWarning("Cannot find controller: meshSet");
        return;
    }

    const auto &meshCont = cont.subController("meshSet");
    if (meshCont.found("regions")) {
        std::string rnl = meshCont.findText("regions");
        replaceAllMarks(rnl, "\n", " ");
        if (!rnl.empty()) {
            size_t pos = 0;
            stdStringList rll;
            split(rnl, rll, ",");
            integer j = 0;
            for (integer i = 0; i < rll.size(); ++i) {
                trim(rll[i]);
                if (!meshCont.found(rll[i])) {
                    checkWarningStr(("Cannot find region: " + rll[i]));
                    continue;
                }
                const auto &regionCont = meshCont.subController(rll[i]);
                regions_.append(markRegion::creator(regionCont));
                regionMap_.emplace(rll[i], regions_(j));
                j++;
            }
        }
    }
}

void OpenHurricane::runtimeMesh::clearRegion() noexcept {
    regions_.clear();
    regionMap_.clear();
}
