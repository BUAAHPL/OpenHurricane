/*!
 * \file rareSlipWall.cpp
 * \brief Main subroutines for rareSlipWall.
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
#include "rareSlipWall.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::rareSlipWall::className_ = "rareSlipWall";

namespace OpenHurricane {
    registerObjFty(vectorBoundary, rareSlipWall, controller);
} // namespace OpenHurricane

OpenHurricane::rareSlipWall::rareSlipWall(const faceZone &fZ, geometryArray<vector, cellMesh> &gf,
                                            const controller &cont)
    : Base(fZ, gf, cont), A_(1), sigma_(1) {
    Base::setSpecified();
    if (cont.found("rareSlip")) {
        A_ = cont.subController("rareSlip").findOrDefault<real>("A", A_);
        sigma_ = cont.subController("rareSlip").findOrDefault<real>("sigma", sigma_);
    }
}

void OpenHurricane::rareSlipWall::updateBoundary() {
    auto &v = this->varField();

    const registerTable &tb = v.tb();

    const auto &T = tb.findObjectRef<cellRealArray>("T");
    const auto &rho = tb.findObjectRef<cellRealArray>("rho");
    const auto &g = tb.findObject<cellRealArray>("gamma");
    const auto &mu = tb.findObject<cellRealArray>("mu");

    const thermoList &thTable = tb.findObject<thermoList>("thermoList");
    const faceList &fL = T.mesh().faces();
    const auto &fA = T.mesh().faceArea();
    const auto &fC = T.mesh().faceCentre();
    const auto &cC = T.mesh().cellCentre();

    const auto &fw = T.mesh().faceWgt();
    realArray yi(thTable.species().size(), Zero);

    for (integer fi = boundaryZone_.firstIndex(); fi < boundaryZone_.lastIndex() + 1; fi++) {
        const integer cl = fL[fi].leftCell();
        const integer cr = fL[fi].rightCell();

        // [m/s]
        const auto &vl = v[cl];
        vector normal = fA[fi].normalized();
        // [m/s]
        real vn = vl * normal;

        // [m]
        const auto d = dist(cC[cl], fC[fi]);

        // [m/s]
        real vt = sqrt(vl * vl - vn * vn);
        // [m/s]
        auto vvt = vl - vn * normal;

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

        const real theta = 2 * A_ * (2 - sigma_) / sigma_ * lambdaw / d;

        const real vss = theta * vt / (1 + theta);

        auto boundaryValue_ = vss * vvt / vt;
        v[cr] = 2.0 * boundaryValue_ - v[cl];
    }
}
