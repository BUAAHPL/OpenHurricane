/*!
 * \file pressureInlet.cpp
 * \brief Main subroutines for pressureInlet.
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
#include "pressureInlet.hpp"
#include "parsingDirection.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::pressureInlet::className_ = "pressureInlet";

namespace OpenHurricane {
    registerObjFty(realBoundary, pressureInlet, controller);
} // namespace OpenHurricane

OpenHurricane::pressureInlet::pressureInlet(const faceZone &fZ, geometryArray<real, cellMesh> &gf,
                                              const controller &cont)
    : Base(fZ, gf, cont), directionType_(normalToBoundary) {
    Base::setSpecified();
    readValue(p0_, "totalPressure", cont);
    readValue(p_, "staticPressure", cont);
    readValue(T0_, "totalTemperature", cont);
    readValue(T_, "T", cont);
    readValue(rho_, "rho", cont);
    readValue(v_, "v", cont);

    const auto dtw = parsingDirection::getDirection(const_cast<controller &>(cont),
                                                    varArray_.mesh(), fZ, direct_);
    if (dtw != "normalToBoundary") {
        directionType_ = givenDirect;
    }
}

void OpenHurricane::pressureInlet::updateBoundary() {
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
        const integer cl = fL[fi].leftCell();
        const integer cr = fL[fi].rightCell();
        vector normal = fA[fi].normalized();
        const real vt = v[cl] * normal;
        const real cd = sqrt(g[cl] * p[cl] / rho[cl]);
        real mad = vt / cd;
        real ma = v[cl].magnitude() / cd;
        real pbnd = p[cl];
        real Tbnd;
        real rhobnd;
        vector vbnd;
        if (mad >= real(1.0)) {
            pbnd = p_;
            Tbnd = T_;
            rhobnd = rho_;
            if (directionType_ == normalToBoundary) {
                vbnd = v_ * normal;
            } else {
                vbnd = v_ * direct_;
            }
        } else {
            pbnd = p0_ / pow(1.0 + (g[cl] - 1.0) / 2.0 * sqr(ma), g[cl] / (g[cl] - 1.0));
            realArray yi(thTable.species().size(), Zero);
            if (thTable.species().size() == 1) {
                yi = real(1.0);
            } else {
                for (integer i = 0; i < thTable.species().size(); i++) {
                    const string &name = thTable.species().name(i);
                    const cellRealArray &sp = tb.findObject<cellRealArray>(name);
                    yi[i] = sp[cl];
                }
            }
            Tbnd = T0_ / (1.0 + (g[cl] - 1.0) / 2.0 * sqr(ma));
            rhobnd = thTable.rho(pbnd, Tbnd, yi);
            if (directionType_ == normalToBoundary) {
                vbnd = ma * cd * normal;
            } else {
                vbnd = ma * cd * direct_;
            }
        }

        if (rhobnd <= 0 || (2.0 * rhobnd - rho[cl] < 0)) {
            pbnd = p[cl];
            Tbnd = T[cl];
            rhobnd = rho[cl];
            vbnd = vt;
        }

        p[cr] = real(2.0) * pbnd - p[cl];
        if (p[cr] < 0) {
            p[cr] = p[cl];
            pbnd = p[cl];
            Tbnd = T[cl];
            rhobnd = rho[cl];
            vbnd = vt;
        }

        rho[cr] = 2.0 * rhobnd - rho[cl];

        v[cr] = real(2.0) * vbnd - v[cl];
        T[cr] = real(2.0) * Tbnd - T[cl];
    }
}
