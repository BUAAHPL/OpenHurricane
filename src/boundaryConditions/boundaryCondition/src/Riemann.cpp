/*!
 * \file Riemann.cpp
 * \brief Main subroutines for Riemann.
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
#include "Riemann.hpp"
#include "parsingDirection.hpp"
#include "registerTable.hpp"
#include "thermoList.hpp"

const std::string OpenHurricane::Riemann::className_ = "Riemann";

namespace OpenHurricane {
    registerObjFty(realBoundary, Riemann, controller);
} // namespace OpenHurricane

OpenHurricane::Riemann::Riemann(const faceZone &fZ, geometryArray<real, cellMesh> &gf,
                                  const controller &cont)
    : Base(fZ, gf, cont) {
    rhoFree_ = cont.findType<real>("rho", rhoFree_);
    pFree_ = cont.findType<real>("p", pFree_);
    tFree_ = cont.findType<real>("T", tFree_);
    gFree_ = cont.findType<real>("gamma", gFree_);

    this->readVector(vFree_, "v", cont);

    cFree_ = sqrt(gFree_ * pFree_ / rhoFree_);
    sFree_ = pFree_ / pow(rhoFree_, gFree_);

    const registerTable &tb = this->varField().tb();
    const thermoList &thTable = tb.findObject<thermoList>("thermoList");
    yiFree_.resize(thTable.species().size(), Zero);
    if (thTable.species().size() == 1) {
        yiFree_ = real(1.0);
    } else {
        for (integer i = 0; i < thTable.species().size(); i++) {
            const string &name = thTable.species().name(i);
            yiFree_[i] = cont.findType<real>(name, yiFree_[i]);
        }
    }
}

void OpenHurricane::Riemann::updateBoundary() {
    const runtimeMesh &mesh = this->mesh();
    const faceZone &fZ = this->boundaryfZ();
    const faceList &fL = mesh.faces();
    const registerTable &tb = this->varField().tb();
    cellRealArray &p = tb.findObjectRef<cellRealArray>("p");
    cellRealArray &rho = tb.findObjectRef<cellRealArray>("rho");
    cellVectorArray &v = tb.findObjectRef<cellVectorArray>("v");
    cellRealArray &T = tb.findObjectRef<cellRealArray>("T");
    const cellRealArray &g = tb.findObject<cellRealArray>("gamma");
    const vectorArray &fa = mesh.faceArea();

    vector vfree = vFree_;
    real co = cFree_;
    const thermoList &thTable = tb.findObject<thermoList>("thermoList");
    realArray yi(thTable.species().size(), Zero);

    for (integer i = fZ.firstIndex(); i < fZ.lastIndex() + 1; i++) {
        integer internalId = fL[i].leftCell();
        integer ghostId = fL[i].rightCell();

        vector fn = -fa[i].normalized();
        real ci = sqrt(g[internalId] * p[internalId] / rho[internalId]);
        real vi = v[internalId] * fn;
        real mi = vi / ci;
        real vo = vfree * fn;
        real g1 = g[internalId] - real(1);

        real rPlus = vi + real(2.0) * ci / g1;
        real rMinus = vo - real(2.0) * co / g1;

        if (mi > real(1.0)) {
            rMinus = vi - real(2.0) * ci / g1;
        } else if (mi < real(-1.0)) {
            rPlus = vo + real(2.0) * co / g1;
        }

        real vb = real(0.5) * (rPlus + rMinus);
        real cb = real(0.25) * g1 * (rPlus - rMinus);

        if (vb < real(0.0)) {
            v[ghostId] = vfree + fn * (vb - vo);

            real sb = sFree_;
            rho[ghostId] = pow(cb * cb / (g[internalId] * sb), real(1.0) / g1);
            p[ghostId] = sb * pow(rho[ghostId], g[internalId]);
            if (thTable.species().size() == 1) {
                yi = real(1.0);
            } else {
                for (integer i = 0; i < thTable.species().size(); i++) {
                    const string &name = thTable.species().name(i);
                    cellRealArray &sp = tb.findObjectRef<cellRealArray>(name);
                    sp[ghostId] = yiFree_[i];
                    yi[i] = sp[ghostId];
                }
            }
        } else {
            v[ghostId] = v[internalId] + fn * (vb - vi);

            real sb = p[internalId] / pow(rho[internalId], g[internalId]);
            rho[ghostId] = pow(cb * cb / (g[internalId] * sb), real(1.0) / g1);
            p[ghostId] = sb * pow(rho[ghostId], g[internalId]);
            if (thTable.species().size() == 1) {
                yi = real(1.0);
            } else {
                for (integer i = 0; i < thTable.species().size(); i++) {
                    const string &name = thTable.species().name(i);
                    cellRealArray &sp = tb.findObjectRef<cellRealArray>(name);
                    sp[ghostId] = sp[internalId];
                    yi[i] = sp[ghostId];
                }
            }
        }
        T[ghostId] = thTable.eos().Tm(rho[ghostId], p[ghostId], yi);
    }
}
