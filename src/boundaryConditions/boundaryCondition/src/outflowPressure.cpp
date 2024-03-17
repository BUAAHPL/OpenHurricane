/*!
 * \file outflowPressure.cpp
 * \brief Main subroutines for outflowPressure.
 * \author Yang Hongzhen
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
#include "outflowPressure.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::outflowPressure::className_ = "outflowPressure";

namespace OpenHurricane {
    registerObjFty(realBoundary, outflowPressure, controller);
} // namespace OpenHurricane

OpenHurricane::outflowPressure::outflowPressure(const faceZone &fZ,
                                                  geometryArray<real, cellMesh> &gf,
                                                  const controller &cont)
    : Base(fZ, gf, cont) {
    this->readValue(p_, std::string("p"), cont);
    this->setSpecified();
}
void OpenHurricane::outflowPressure::updateBoundary() {
    const registerTable &tb = this->varField().tb();
    cellRealArray &p = this->varField();
    cellVectorArray &v = tb.findObjectRef<cellVectorArray>("v");
    cellRealArray &rho = tb.findObjectRef<cellRealArray>("rho");
    cellRealArray &T = tb.findObjectRef<cellRealArray>("T");
    const cellRealArray &g = tb.findObject<cellRealArray>("gamma");

    const thermoList &thTable = tb.findObject<thermoList>("thermoList");
    const faceList &fL = p.mesh().faces();
    const vectorArray &fA = p.mesh().faceArea();

    for (integer fi = boundaryZone_.firstIndex(); fi < boundaryZone_.lastIndex() + 1; fi++) {
        integer celli = fL[fi].leftCell();
        vector normal = -fA[fi].normalized();

        real ma = v[celli] * normal / sqrt(g[celli] * p[celli] / rho[celli]);

        integer ghostId = fL[fi].rightCell();
        if (ma >= real(1.0)) {
            p[ghostId] = p[celli];

        } else {
            p[ghostId] = p_;
        }
    }
}