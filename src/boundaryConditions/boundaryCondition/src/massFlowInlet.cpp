/*!
 * \file massFlowInlet.cpp
 * \brief Main subroutines for massFlowInlet.
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
#include "massFlowInlet.hpp"
#include "parsingDirection.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::massFlowInlet::className_ = "massFlowInlet";

namespace OpenHurricane {
    registerObjFty(realBoundary, massFlowInlet, controller);
} // namespace OpenHurricane

OpenHurricane::massFlowInlet::massFlowInlet(const faceZone &fZ, geometryArray<real, cellMesh> &gf,
                                              const controller &cont)
    : Base(fZ, gf, cont), directionType_(normalToBoundary) {
    Base::setSpecified();
    this->readValue(h0_, std::string("totalEnthalpy"), cont);
    this->readValue(Tt_, std::string("totalTemperature"), cont);
    this->readValue(flux_, std::string("massFlux"), cont);
    this->readValue(p_, std::string("p"), cont);
    const auto dtw = parsingDirection::getDirection(const_cast<controller &>(cont),
                                                    varArray_.mesh(), fZ, direct_);
    if (dtw != "normalToBoundary") {
        directionType_ = givenDirect;
    }

    const registerTable &tb = this->varField().tb();
    const thermoList &thTable = tb.findObject<thermoList>("thermoList");
    yiBnd_.resize(thTable.species().size(), real(1));
    for (integer i = 0; i < thTable.species().size(); ++i) {
        std::string spn = thTable.species()[i].name();
        yiBnd_[i] = cont.findOrDefault<real>(spn, 0.0);
    }
}

void OpenHurricane::massFlowInlet::updateBoundary() {
    const registerTable &tb = this->varField().tb();
    cellRealArray &rho = this->varField();
    cellVectorArray &v = tb.findObjectRef<cellVectorArray>("v");
    cellRealArray &p = tb.findObjectRef<cellRealArray>("p");
    cellRealArray &T = tb.findObjectRef<cellRealArray>("T");
    const cellRealArray &g = tb.findObject<cellRealArray>("gamma");

    const thermoList &thTable = tb.findObject<thermoList>("thermoList");
    const faceList &fL = rho.mesh().faces();
    const vectorArray &fA = rho.mesh().faceArea();

    const auto &fw = rho.mesh().faceWgt();
    real Testd = 0;

    real h0 = thTable.ha_p(p_, Tt_, yiBnd_);

    for (integer fi = boundaryZone_.firstIndex(); fi < boundaryZone_.lastIndex() + 1; fi++) {
        const integer celli = fL[fi].leftCell();
        const integer cr = fL[fi].rightCell();

        vector normal = fA[fi].normalized();
        real vt = v[celli] * normal;
        real ma = vt / sqrt(g[celli] * p[celli] / rho[celli]);
        real pbnd = p[celli];
        if (ma >= real(1.0)) {
            pbnd = p_;
        }

        if (fi == boundaryZone_.firstIndex()) {
            Testd = T[celli];
        }

        real Tnewd = Testd;
        real Ttold = Tnewd * real(1.0e-5);
        int iter = 0;
        real pd = pbnd;
        do {
            Testd = Tnewd;
            real rhobndd = thTable.rho(pd, Testd, yiBnd_);
            real vtd = flux_ / rhobndd;

            real dFd = thTable.ha_p(pd, Testd, yiBnd_) + real(0.5) * vtd * vtd - h0;

            real dFdTd = thTable.cp_p(pd, Testd, yiBnd_) +
                         vtd * flux_ * thTable.eos().DRhoInverseDT(pd, rhobndd, Testd, yiBnd_);
            real dTd = dFd / dFdTd;
            integer flag = 1;
            Tnewd = thTable.limit(Testd - dTd, flag);
            if (iter++ > 100) {
#ifdef HUR_DEBUG
                PLWarning("Maximum number of iterations exceeded: 100");
#else
                LWarning("Maximum number of iterations exceeded: 100");
#endif // HUR_DEBUG
                break;
            }
        } while (mag(Tnewd - Testd) > Ttold);
        real Tbnd = Tnewd;
        const real rhobnd = thTable.rho(pd, Tnewd, yiBnd_);
        bool negativeGhostDensity = false;
        real boundaryValue_ = 0;
        if (2.0 * rhobnd - rho[celli] < 0.0) {
            negativeGhostDensity = true;
            boundaryValue_ = 0.5 * (rhobnd + rho[celli]);
        } else {
            boundaryValue_ = rhobnd;
        }

        rho[cr] = 2.0 * boundaryValue_ - rho[celli];

        const real vnbnd = flux_ / rhobnd;
        vector vbnd;
        if (directionType_ == normalToBoundary) {
            vbnd = vnbnd * normal;
        } else {
            vbnd = vnbnd * direct_;
        }

        integer ghostId = fL[fi].rightCell();

        if (negativeGhostDensity) {
            v[ghostId] = vbnd;
            p[ghostId] = pbnd;
            T[ghostId] = Tbnd;
        } else {
            v[ghostId] = real(2.0) * vbnd - v[celli];
            p[ghostId] = real(2.0) * pbnd - p[celli];
            T[ghostId] = real(2.0) * Tbnd - T[celli];
        }
    }
}
