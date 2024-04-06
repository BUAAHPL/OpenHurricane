/*!
 * \file pressureOutlet.cpp
 * \brief Main subroutines for pressureOutlet.
 * \author Yang Hongzhen
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
#include "pressureOutlet.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::pressureOutlet::className_ = "pressureOutlet";

namespace OpenHurricane {
    registerObjFty(realBoundary, pressureOutlet, controller);
} // namespace OpenHurricane

OpenHurricane::pressureOutlet::pressureOutlet(const faceZone &fZ,
                                                geometryArray<real, cellMesh> &gf,
                                                const controller &cont)
    : Base(fZ, gf, cont) {
    this->readValue(p_, std::string("p"), cont);
    this->setSpecified();
}

void OpenHurricane::pressureOutlet::updateBoundary() {
    const registerTable &tb = this->varField().tb();
    cellRealArray &rho = this->varField();
    cellVectorArray &v = tb.findObjectRef<cellVectorArray>("v");
    cellRealArray &p = tb.findObjectRef<cellRealArray>("p");
    cellRealArray &T = tb.findObjectRef<cellRealArray>("T");
    const cellRealArray &g = tb.findObject<cellRealArray>("gamma");

    const thermoList &thTable = tb.findObject<thermoList>("thermoList");
    const faceList &fL = this->mesh().faces();
    const vectorArray &fA = this->mesh().faceArea();

    for (integer fi = boundaryZone_.firstIndex(); fi < boundaryZone_.lastIndex() + 1; fi++) {
        integer celli = fL[fi].leftCell();
        vector normal = -fA[fi].normalized();

        real ma = v[celli] * normal / sqrt(g[celli] * p[celli] / rho[celli]);
        real Tbnd = T[celli];
        realArray yi(thTable.species().size(), Zero);
        if (thTable.species().size() == 1) {
            yi = real(1.0);
        } else {
            for (integer i = 0; i < thTable.species().size(); i++) {
                const string &name = thTable.species().name(i);
                const cellRealArray &sp = tb.findObject<cellRealArray>(name);
                yi[i] = sp[celli];
            }
        }
        real pbnd = p_;
        if (ma > real(1.0)) {
            pbnd = p[celli];
        }

        real rhobnd = thTable.rho(pbnd, Tbnd, yi);

        integer ghostId = fL[fi].rightCell();
        p[ghostId] = real(2.0) * pbnd - p[celli];
        if (p[ghostId] < 0) {
            pbnd = p[celli];
            p[ghostId] = p[celli];
            Tbnd = T[celli];
            rhobnd = rho[celli];
        }
        rho[ghostId] = 2.0 * rhobnd - rho[celli];
    }
}