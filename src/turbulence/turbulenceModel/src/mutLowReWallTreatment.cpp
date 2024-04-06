/*!
 * \file mutLowReWallTreatment.cpp
 * \brief Main subroutines of mut near-wall low-Re treatment.
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
#include "mutLowReWallTreatment.hpp"

namespace OpenHurricane {
    createClassNameStr(mutLowReWallTreatment, "mutLowReWallTreatment");
    registerObjFty(realBoundary, mutLowReWallTreatment, controller);
} // namespace OpenHurricane

OpenHurricane::mutLowReWallTreatment::mutLowReWallTreatment(const faceZone &fZ,
                                                            geometryArray<real, cellMesh> &gf,
                                                            const controller &cont)
    : Base(fZ, gf, cont), mutBnd_(cont.findOrDefault<real>("mut", 0)) {}

OpenHurricane::mutLowReWallTreatment::mutLowReWallTreatment(const mutLowReWallTreatment &bB)
    : Base(bB), mutBnd_(bB.mutBnd_) {}

void OpenHurricane::mutLowReWallTreatment::updateBoundary() {
    const registerTable &tb = this->varField().tb();

    const auto &faces = mesh().faces();
    auto &mut = varArray_;
    for (integer facei = boundaryZone_.firstIndex(); facei <= boundaryZone_.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        const integer cr = faces[facei].rightCell();
        mut[cr] = 2.0 * mutBnd_ - mut[cl];
    }
}