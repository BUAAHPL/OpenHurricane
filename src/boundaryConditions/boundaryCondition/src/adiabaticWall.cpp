/*!
 * \file adiabaticWall.cpp
 * \brief Main subroutines for adiabaticWall.
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
#include "adiabaticWall.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::adiabaticWall::className_ = "adiabaticWall";

namespace OpenHurricane {
    registerObjFty(realBoundary, adiabaticWall, controller);
} // namespace OpenHurricane

OpenHurricane::adiabaticWall::adiabaticWall(const faceZone &fZ, geometryArray<real, cellMesh> &gf,
                                              const controller &cont)
    : Base(fZ, gf, cont) {
    Base::setSpecified();
    const registerTable &tb = this->varField().tb();
    const thermoList &thTable = tb.findObject<thermoList>("thermoList");

    yiBnd_.resize(thTable.species().size(), Zero);
}

void OpenHurricane::adiabaticWall::updateBoundary() {
    const registerTable &tb = this->varField().tb();
    cellRealArray &rho = tb.findObjectRef<cellRealArray>("rho");
    cellRealArray &p = tb.findObjectRef<cellRealArray>("p");
    cellRealArray &T = this->varField();

    const thermoList &thTable = tb.findObject<thermoList>("thermoList");

    const faceList &fL = T.mesh().faces();
    const auto &fw = T.mesh().faceWgt();

    auto &yi = yiBnd_;
    for (integer fi = boundaryZone_.firstIndex(); fi < boundaryZone_.lastIndex() + 1; fi++) {       
        const auto cl = fL[fi].leftCell();
        const auto cr = fL[fi].rightCell();
        T[cr] = T[cl];

        if (thTable.species().size() == 1) {
            yi = real(1.0);
        } else {
            for (integer i = 0; i < thTable.species().size(); i++) {
                const string &name = thTable.species().name(i);
                const cellRealArray &sp = tb.findObject<cellRealArray>(name);
                yi[i] = fw[fi] * sp[cl] + (1 - fw[fi]) * sp[cr];
            }
        }

        rho[cr] = thTable.rho(p[cl], T[cr], yi);
    }
}
