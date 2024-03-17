/*!
 * \file mutKOmegaInlet.cpp
 * \brief Main subroutines of mut k-omega inlet boundary.
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
#include "mutKOmegaInlet.hpp"

namespace OpenHurricane {
    createClassNameStr(mutKOmegaInlet, "mutKOmegaInlet");
    registerObjFty(realBoundary, mutKOmegaInlet, controller);
} // namespace OpenHurricane

OpenHurricane::mutKOmegaInlet::mutKOmegaInlet(const faceZone &fZ, geometryArray<real, cellMesh> &gf,
                                              const controller &cont)
    : Base(fZ, gf, cont), cv1_(7.1) {}

OpenHurricane::mutKOmegaInlet::mutKOmegaInlet(const mutKOmegaInlet &bB) : Base(bB), cv1_(bB.cv1_) {}

void OpenHurricane::mutKOmegaInlet::updateBoundary() {
    const registerTable &tb = this->varField().tb();
    const auto &rho = tb.findObjectRef<cellRealArray>("rho");
    const auto &kt = tb.findObjectRef<cellRealArray>("kt");
    const auto &wt = tb.findObjectRef<cellRealArray>("wt");

    const runtimeMesh &mesh = this->mesh();
    const faceZone &fZ = this->boundaryfZ();
    const faceList &fL = mesh.faces();
    auto &mut = this->varArray_;

    for (integer i = fZ.firstIndex(); i < fZ.lastIndex() + 1; i++) {
        const integer internalId = fL[i].leftCell();
        const integer ghostId = fL[i].rightCell();
        mut[ghostId] = rho[ghostId] * kt[ghostId] / wt[ghostId];
    }
}
