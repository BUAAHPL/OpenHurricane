/*!
 * \file omegaNearWallTreatment.cpp
 * \brief Main subroutines of omega near-wall treatment.
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
#include "nutWallTreatment.hpp"

namespace OpenHurricane {
    createClassNameStr(nutWallTreatment, "nutWallTreatment");
    registerObjFty(realBoundary, nutWallTreatment, controller);
} // namespace OpenHurricane

OpenHurricane::nutWallTreatment::nutWallTreatment(const faceZone &fZ,
                                                  geometryArray<real, cellMesh> &gf,
                                                  const controller &cont)
    : Base(fZ, gf, cont), nutBnd_(cont.findOrDefault<real>("nut", 0)) {}

OpenHurricane::nutWallTreatment::nutWallTreatment(const nutWallTreatment &bB)
    : Base(bB), nutBnd_(bB.nutBnd_) {}

void OpenHurricane::nutWallTreatment::updateBoundary() {
    const registerTable &tb = this->varField().tb();

    const auto &faces = mesh().faces();
    auto &nut = varArray_;
    for (integer facei = boundaryZone_.firstIndex(); facei <= boundaryZone_.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        const integer cr = faces[facei].rightCell();
        nut[cr] = 2.0 * nutBnd_ - nut[cl];
    }
}
