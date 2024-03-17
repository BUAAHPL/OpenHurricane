/*!
 * \file cellBasedMethod.cpp
 * \brief Main subroutines for cell-based method for computing time step.
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

#include "cellBasedMethod.hpp"
#include "cellArrays.hpp"

namespace OpenHurricane {
    createClassNameStr(cellBasedMethod, "cellBasedMethod");
    registerObjFty(pseudoTime, cellBasedMethod, controller);
} // namespace OpenHurricane

OpenHurricane::realArray &OpenHurricane::cellBasedMethod::deltaSx() {
    if (!deltaSxPtr_) {
        deltaSxPtr_.reset(new realArray(mesh_.nCells()));
        const auto &faces = mesh_.faces();
        const auto &cells = mesh_.cells();
        const auto &fA = mesh_.faceArea();
        for (integer n = 0; n < mesh_.nCells(); ++n) {
            real ra_sum = 0.0;
            (*deltaSxPtr_)[n] = Zero;
            for (integer i = 0; i < cells[n].faceSize(); ++i) {
                const integer j = cells[n].facei(i);
                const auto &cl = faces[j].leftCell();
                const auto &cr = faces[j].rightCell();
                (*deltaSxPtr_)[n] += mag(fA[j].x());
            }
            (*deltaSxPtr_)[n] *= 0.5;
        }
    }
    return *deltaSxPtr_;
}

OpenHurricane::realArray &OpenHurricane::cellBasedMethod::deltaSy() {
    if (!deltaSyPtr_) {
        deltaSyPtr_.reset(new realArray(mesh_.nCells()));
        const auto &faces = mesh_.faces();
        const auto &cells = mesh_.cells();
        const auto &fA = mesh_.faceArea();
        for (integer n = 0; n < mesh_.nCells(); ++n) {
            real ra_sum = 0.0;
            (*deltaSyPtr_)[n] = Zero;
            for (integer i = 0; i < cells[n].faceSize(); ++i) {
                const integer j = cells[n].facei(i);
                const auto &cl = faces[j].leftCell();
                const auto &cr = faces[j].rightCell();
                (*deltaSyPtr_)[n] += mag(fA[j].y());
            }
            (*deltaSyPtr_)[n] *= 0.5;
        }
    }
    return *deltaSyPtr_;
}

OpenHurricane::realArray &OpenHurricane::cellBasedMethod::deltaSz() {
    if (!deltaSzPtr_) {
        deltaSzPtr_.reset(new realArray(mesh_.nCells()));
        const auto &faces = mesh_.faces();
        const auto &cells = mesh_.cells();
        const auto &fA = mesh_.faceArea();
        for (integer n = 0; n < mesh_.nCells(); ++n) {
            real ra_sum = 0.0;
            (*deltaSzPtr_)[n] = Zero;
            for (integer i = 0; i < cells[n].faceSize(); ++i) {
                const integer j = cells[n].facei(i);
                const auto &cl = faces[j].leftCell();
                const auto &cr = faces[j].rightCell();
                (*deltaSzPtr_)[n] += mag(fA[j].z());
            }
            (*deltaSzPtr_)[n] *= 0.5;
        }
    }
    return *deltaSzPtr_;
}

OpenHurricane::cellBasedMethod::cellBasedMethod(
    const controller &cont, const runtimeMesh &mesh, const iteration &iter, const flowModel &flows,
    realArray &dt, const cellRealArray &shockFactor, const faceVector2DArray &rai,
    const faceVector2DArray &rav, const integerArray &temperatureFlag,
    const integerArray &pressureFlag, const cellIntegerArray &CFLFlag)
    : pseudoTime(cont, mesh, iter, flows, dt, shockFactor, rai, rav, temperatureFlag, pressureFlag,
                 CFLFlag),
      deltaSxPtr_(nullptr), deltaSyPtr_(nullptr), deltaSzPtr_(nullptr), useLowMachPrecon_(false) {}

OpenHurricane::cellBasedMethod::~cellBasedMethod() noexcept {
    deltaSxPtr_.clear();
    deltaSyPtr_.clear();
    deltaSzPtr_.clear();
}

void OpenHurricane::cellBasedMethod::computingTimeStep() {
    computingTimeStep(dt_, cflPtr_->getCFL());
}

void OpenHurricane::cellBasedMethod::computingTimeStep(realArray &dt, const real cfl0) {

    if (mesh_.moving()) {
        deltaSxPtr_.clear();
        deltaSyPtr_.clear();
        deltaSzPtr_.clear();
    }
    const real C = CForTimeStep_;

    const auto &faces = mesh_.faces();
    const auto &cells = mesh_.cells();
    const auto &cV = mesh_.cellVolume();
    const auto &fA = mesh_.faceArea();

    const auto &v = flow_.v();
    const auto &rho = flow_.rho();
    const auto &p = flow_.p();

    const auto &T = flow_.T();
    const auto &mu = flow_.mul();
    const auto &g = flow_.gama();
    const auto &dSx = deltaSx();
    const auto &dSy = deltaSy();
    const auto &dSz = deltaSz();
    if (useLowMachPrecon_) {
        const auto &a4 = flow_.mesh().findObjectRef<cellRealArray>("a4ForLowMachPrecond");
        const auto &a5 = flow_.mesh().findObjectRef<cellRealArray>("a5ForLowMachPrecond");
        for (integer n = 0; n < mesh_.nCells(); ++n) {
            const real vv = v[n].magnitude();
            real ffa = 0.5 * (a4[n] + 1.0);
            real c = 0.5 * sqrt(sqr(vv) * sqr(a4[n] - 1.0) + real(4) * a5[n]);

            const real Acx = (ffa * mag(v[n].x()) + c) * dSx[n];
            const real Acy = (ffa * mag(v[n].y()) + c) * dSy[n];
            const real Acz = (ffa * mag(v[n].z()) + c) * dSz[n];

            const real mid = max(real(4.0 / 3.0), g[n]) * flow_.keEff(n) / rho[n] / cV[n];
            const real Avx = mid * sqr(dSx[n]);
            const real Avy = mid * sqr(dSy[n]);
            const real Avz = mid * sqr(dSz[n]);

            dt[n] = cfl0 * cV[n] / (Acx + Acy + Acz + C * (Avx + Avy + Avz));
        }
    } else {
        for (integer n = 0; n < mesh_.nCells(); ++n) {
            const real c = sqrt(g[n] * p[n] / rho[n]);
            const real Acx = (mag(v[n].x()) + c) * dSx[n];
            const real Acy = (mag(v[n].y()) + c) * dSy[n];
            const real Acz = (mag(v[n].z()) + c) * dSz[n];

            const real mid = max(real(4.0 / 3.0), g[n]) * flow_.keEff(n) / rho[n] / cV[n];
            const real Avx = mid * sqr(dSx[n]);
            const real Avy = mid * sqr(dSy[n]);
            const real Avz = mid * sqr(dSz[n]);

            dt[n] = cfl0 * cV[n] / (Acx + Acy + Acz + C * (Avx + Avy + Avz));
        }
    }
}
