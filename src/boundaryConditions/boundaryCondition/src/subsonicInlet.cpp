/*!
 * \file subsonicInlet.cpp
 * \brief Main subroutines for subsonicInlet.
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
#include "subsonicInlet.hpp"
#include "parsingDirection.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::subsonicInlet::className_ = "subsonicInlet";

namespace OpenHurricane {
    registerObjFty(realBoundary, subsonicInlet, controller);
} // namespace OpenHurricane

OpenHurricane::subsonicInlet::subsonicInlet(const faceZone &fZ, geometryArray<real, cellMesh> &gf,
                                              const controller &cont)
    : Base(fZ, gf, cont), directionType_(normalToBoundary) {
    Base::setSpecified();
    readValue(p0_, "totalPressure", cont);
    readValue(p_, "staticPressure", cont);
    readValue(T0_, "totalTemperature", cont);

    this->setSpecified();
    const auto dtw =
        parsingDirection::getDirection(const_cast<controller &>(cont), mesh(), fZ, direct_);
    if (dtw != "normalToBoundary") {
        directionType_ = givenDirect;
    }
}

void OpenHurricane::subsonicInlet::updateBoundary() {
    const registerTable &tb = this->varField().tb();
    cellRealArray &rho = this->varField();
    cellVectorArray &v = tb.findObjectRef<cellVectorArray>("v");
    cellRealArray &p = tb.findObjectRef<cellRealArray>("p");
    cellRealArray &T = tb.findObjectRef<cellRealArray>("T");
    const cellRealArray &g = tb.findObject<cellRealArray>("gamma");

    const thermoList &thTable = tb.findObject<thermoList>("thermoList");
    const faceList &fL = this->mesh().faces();
    const auto &fA = this->mesh().faceArea();

    for (integer fi = boundaryZone_.firstIndex(); fi < boundaryZone_.lastIndex() + 1; fi++) {
        const integer cl = fL[fi].leftCell();
        const integer cr = fL[fi].rightCell();

        const real gamma = g[cl];

        const real cd = sqrt(gamma * p[cl] / rho[cl]);
        const vector n = -fA[fi].normalized();
        const real vdn = v[cl] * n;
        const real vd = v[cl].magnitude();
        const real rm = vdn - 2.0 * cd / (gamma - 1.0);
        const real cosTheta = -vdn / vd;
        const real c02 = sqr(cd) + (gamma - 1.0) / 2.0 * sqr(vd);
        const real cb = -rm * (gamma - 1.0) / ((gamma - 1.0) * sqr(cosTheta) + 2.0) *
                        (1.0 + cosTheta * sqrt(((gamma - 1.0) * sqr(cosTheta) + 2.0) * c02 /
                                                   ((gamma - 1.0) * sqr(rm)) -
                                               (gamma - 1.0) / 2.0));
        const real Tbnd = T0_ * (sqr(cb) / c02);
        const real pbnd = p0_ * pow((Tbnd / T0_), gamma / (gamma - real(1.0)));
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
        const real rhobnd = thTable.rho(pbnd, Tbnd, yi);
        const real cpbnd = thTable.cp_p(pbnd, Tbnd, yi);
        const real vbndn = sqrt(2.0 * cpbnd * (T0_ - Tbnd));
        vector vbnd;
        if (directionType_ == normalToBoundary) {
            vbnd = -vbndn * n;
        } else {
            vbnd = vbndn * direct_;
        }
        rho[cr] = real(2.0) * rhobnd - rho[cl];

        v[cr] = real(2.0) * vbnd - v[cl];
        p[cr] = real(2.0) * pbnd - p[cl];
        T[cr] = real(2.0) * Tbnd - T[cl];
    }
}
