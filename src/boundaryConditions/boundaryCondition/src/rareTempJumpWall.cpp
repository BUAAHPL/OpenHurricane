/*!
 * \file rareTempJumpWall.cpp
 * \brief Main subroutines for rareTempJumpWall.
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
#include "rareTempJumpWall.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::rareTempJumpWall::className_ = "rareTempJumpWall";

namespace OpenHurricane {
    registerObjFty(realBoundary, rareTempJumpWall, controller);
} // namespace OpenHurricane

OpenHurricane::rareTempJumpWall::rareTempJumpWall(const faceZone &fZ,
                                                    geometryArray<real, cellMesh> &gf,
                                                    const controller &cont)
    : Base(fZ, gf, cont), alpha_(1), Twd_(300), Pr_(0.72) {
    Base::setSpecified();
    alpha_ = cont.findOrDefault<real>("alpha", alpha_);
    readValue(Twd_, "T", cont);

    if (cont.topCont().found("flow")) {
        Pr_ = cont.topCont().subController("flow").findOrDefault<real>("Prl", Pr_);
    }
}

void OpenHurricane::rareTempJumpWall::updateBoundary() {
    auto &T = this->varField();

    const registerTable &tb = T.tb();

    auto &rho = tb.findObjectRef<cellRealArray>("rho");
    const auto &g = tb.findObject<cellRealArray>("gamma");
    const auto &mu = tb.findObject<cellRealArray>("mu");
    const cellRealArray &p = tb.findObject<cellRealArray>("p");

    const thermoList &thTable = tb.findObject<thermoList>("thermoList");
    const faceList &fL = T.mesh().faces();
    const auto &fC = T.mesh().faceCentre();
    const auto &cC = T.mesh().cellCentre();

    const auto &fw = T.mesh().faceWgt();
    realArray yi(thTable.species().size(), Zero);

    for (integer fi = boundaryZone_.firstIndex(); fi < boundaryZone_.lastIndex() + 1; fi++) {
        const integer cl = fL[fi].leftCell();
        const integer cr = fL[fi].rightCell();

        // [m/s]
        const auto &Tl = T[cl];

        // [m]
        const auto d = dist(cC[cl], fC[fi]);

        // [kg/(m s)]
        const real muw = (fw[fi] * mu[cl] + (1 - fw[fi]) * mu[cr]);
        const real gammaw = fw[fi] * g[cl] + (1 - fw[fi]) * g[cr];

        // [kg/m^3]
        const real rhow = (fw[fi] * rho[cl] + (1 - fw[fi]) * rho[cr]);

        // [K]
        const real Tw = (fw[fi] * T[cl] + (1 - fw[fi]) * T[cr]);
        if (thTable.species().size() == 1) {
            yi = real(1.0);
        } else {
            for (integer i = 0; i < thTable.species().size(); i++) {
                const string &name = thTable.species().name(i);
                const cellRealArray &sp = tb.findObject<cellRealArray>(name);
                yi[i] = sp[cl] * fw[fi] + sp[cr] * (1.0 - fw[fi]);
            }
        }
        // [j/kg/K]
        const auto Rw = thTable.species().Rm(yi);

        // [m]
        const real lambdaw = muw / rhow * sqrt(constant::mathConstants::pi / (real(2) * Rw * Tw));

        const real elsilon =
            (2 - alpha_) / alpha_ * 2 * gammaw / (gammaw + 1) / Pr_ * 2 * lambdaw / d;

        const real T0 = (Twd_ + elsilon * Tl) / (1 + elsilon);

        T[cr] = 2.0 * T0 - T[cl];

        auto rhoBnd = thTable.rho(p[cl], T0, yi);

        rho[cr] = 2.0 * rhoBnd - rho[cl];
    }
}
