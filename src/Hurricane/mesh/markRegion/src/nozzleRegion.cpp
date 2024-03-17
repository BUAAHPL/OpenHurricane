/*!
 * \file nozzleRegion.cpp
 * \brief Main subroutines for nozzle region of mesh.
 * \author Peng Jian
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

#include "nozzleRegion.hpp"
#include "cellMesh.hpp"
#include "formulaParsing.hpp"
#include "constants.hpp"

namespace OpenHurricane {
    createClassNameStr(nozzleRegion, "nozzle");
}
namespace OpenHurricane {
    registerObjFty(markRegion, nozzleRegion, controller);
}

hur_nodiscard OpenHurricane::real OpenHurricane::nozzleRegion::divePresFromDia(const real diaRatio) const {
    real squareD = sqr(diaRatio);
    real midValue = 0;
    real testSquareD = 0;
    real minVar = 0;
    real maxVar = criticalValue_;
    real flag = 0;
    for (integer i = 0; i < 100000; i++) {
        midValue = (minVar + maxVar) / 2.0;
        testSquareD = areaRatioFromPres(midValue);
        flag = abs(squareD - testSquareD);
        if (flag < 1e-6) {
            break;
        } else {
            if (squareD > testSquareD) {
                maxVar = midValue;
            } else {
                minVar = midValue;
            }
        }
    }
    return midValue;
}

hur_nodiscard OpenHurricane::real OpenHurricane::nozzleRegion::convPresFromDia(const real diaRatio) const {
    real squareD = sqr(diaRatio);
    real midValue = 0;
    real testSquareD = 0;
    real minVar = criticalValue_;
    real maxVar = real(1.0);
    real flag = 0;
    for (integer i = 0; i < 100000; i++) {
        midValue = (minVar + maxVar) / 2.0;
        testSquareD = areaRatioFromPres(midValue);
        flag = abs(squareD - testSquareD);
        if (flag < 1e-6) {
            break;
        } else {
            if (squareD > testSquareD) {
                minVar = midValue;
            } else {
                maxVar = midValue;
            }
        }
    }
    return midValue;
}

hur_nodiscard OpenHurricane::real
OpenHurricane::nozzleRegion::areaRatioFromPres(const real presRatio) const {
    const real gammaMinus = gammaConst_ - real(1.0);
    const real gammaPlus = gammaConst_ + real(1.0);
    const real m = 2 / gammaPlus;
    const real n = 0.5 * gammaPlus / gammaMinus;
    const real r = sqrt(2 * gammaConst_ / gammaMinus);
    const real f = gammaMinus / gammaConst_;
    const real u = pow(m, n) * sqrt(gammaConst_) / r;
    const real w = pow(presRatio, 1 / gammaConst_) * sqrt(1 - pow(presRatio, f));
    const real areaRatio = u / w;
    return areaRatio;
}

void OpenHurricane::nozzleRegion::calTemperature(realGeometryArray<cellMesh> &T) const {}

void OpenHurricane::nozzleRegion::calDensity(realGeometryArray<cellMesh> &rho) const {}

void OpenHurricane::nozzleRegion::calPressure(realGeometryArray<cellMesh> &p) const {
    if (!pressPatchedFlag_) {
        const auto &mesh = p.mesh();

        const vector Ae = Aexit_;
        const vector Ah = Ahead_;
        const vector At = Athroat_;
        const real De = Dexit_;
        const real Dh = Dhead_;
        const real Dt = Dthroat_;

        const real AtAe = dist(At, Ae);
        const real AtAh = dist(At, Ah);

        const vector throatToHead = ((Ahead_ - Athroat_)).normalized();
        const vector throatToExit = ((Aexit_ - Athroat_)).normalized();

        const auto &cC = mesh.cellCentre();

        for (integer n = 0; n < mesh.nCells(); ++n) {
            const auto &C = cC[n];

            const vector ct = At - C;
            const vector ch = Ah - C;
            const vector ce = Ae - C;
            real presRatio = 0;
            // convergence
            if (ch * throatToHead >= 0 && ct * throatToHead <= 0) {
                vector d = ((At - Ah) ^ (C - Ah)) / AtAh;
                real nLen = d.magnitude();
                real CABA = (C - Ah) * (At - Ah);
                real CBAB = (C - At) * (Ah - At);
                if (isInside()) {
                    if (nLen > Dh)
                        continue;
                    if (CABA < 0.0 || CBAB < 0.0)
                        continue;
                } else {
                    if (nLen <= Dh && CABA > 0.0 && CBAB >= 0.0)
                        continue;
                }
                const real conProject = -ct * throatToHead;
                const real Dc = real(1.0) + conProject * (Dh / Dt - real(1.0)) / AtAh;
                presRatio = convPresFromDia(Dc);
            }

            // divergence
            if (ce * throatToExit >= 0 && ct * throatToExit <= 0) {
                vector d = ((Ae - At) ^ (C - At)) / AtAe;
                real nLen = d.magnitude();
                real CABA = (C - At) * (Ae - At);
                real CBAB = (C - Ae) * (At - Ae);
                if (isInside()) {
                    if (nLen > De)
                        continue;
                    if (CABA < 0.0 || CBAB < 0.0)
                        continue;
                } else {
                    if (nLen <= De && CABA > 0.0 && CBAB >= 0.0)
                        continue;
                }
                const real divProject = -ct * throatToExit;
                const real Dc = real(1.0) + divProject * (De / Dt - real(1.0)) / AtAe;
                presRatio = divePresFromDia(Dc);
            }
            p[n] = presRatio * pc_;
        }
    }
    return;
}

OpenHurricane::nozzleRegion::nozzleRegion(const controller &cont)
    : markRegion(cont), Aexit_(cont.findOrDefault<vector>("Aexit", vector(0.0))),
      Ahead_(cont.findOrDefault<vector>("Ahead", vector(0.0))),
      Athroat_(cont.findOrDefault<vector>("Athroat", vector(0.0))),
      Dhead_(cont.findOrDefault<real>("Dhead", real(0.0))),
      Dexit_(cont.findOrDefault<real>("Dexit", real(0.0))),
      Dthroat_(cont.findOrDefault<real>("Dthroat", real(0.0))), gammaConst_(real(0.0)),
      criticalValue_(real(0.0)), pc_(real(-1.0)), boundaryFaceName_(), pressPatchedFlag_(false),
      varPatchFlagMap_(nullptr) {
    criticalValue_ =
        pow((2 / (gammaConst_ + real(1.0))), (gammaConst_ / (gammaConst_ - real(1.0))));

    if (cont.found("boundaryFaceName")) {
        boundaryFaceName_ = cont.findWord("boundaryFaceName");
    } else {
        LFatal("Boundary face name is essential in nozzle region");
    }
    const vector AtAh = Ahead_ - Athroat_;
    const vector AtAe = Aexit_ - Athroat_;

    if (AtAh * AtAe > 0) {
        LFatal("Invalid nozzle Region settings, please check the "
               "position of those three cross-section");
    }
}

inline void OpenHurricane::nozzleRegion::setGammaConst(const real gammabc) {
    gammaConst_ = gammabc;
}

hur_nodiscard OpenHurricane::integerList
OpenHurricane::nozzleRegion::regionCellId(const runtimeMesh &mesh) const {
    integerList cid;
    if (!isOptionSet()) {
        PLWarning("Option unset for mark region: %d", id());
        return cid;
    }

    vector Ae = Aexit_;
    vector Ah = Ahead_;
    vector At = Athroat_;
    real De = Dexit_;
    real Dh = Dhead_;
    real Dt = Dthroat_;

    real AhAe = dist(Ah, Ae);

    const auto &cC = mesh.cellCentre();
    cid.resize(mesh.nCells(), -1);
    integer count = 0;

    for (integer n = 0; n < mesh.nCells(); ++n) {
        const auto &C = cC[n];

        vector d = ((Ae - Ah) ^ (C - Ah)) / AhAe;
        real nLen = d.magnitude();

        // AC * AB
        real CABA = (C - Ah) * (Ae - Ah);

        // BC * BA
        real CBAB = (C - Ae) * (Ah - Ae);

        if (isInside()) {
            if (nLen > De)
                continue;
            if (CABA < 0.0 || CBAB < 0.0)
                continue;
        } else {
            if (nLen <= De && CABA > 0.0 && CBAB >= 0.0)
                continue;
        }

        cid[count] = n;
        count++;
    }

    cid.resize(count);
    return cid;
}

bool OpenHurricane::nozzleRegion::patching(realGeometryArray<cellMesh> &cellQ, real &value) const {
    return patch(cellQ, value);
}

bool OpenHurricane::nozzleRegion::patching(vectorGeometryArray<cellMesh> &cellQ, vector &value) const {
    return patch(cellQ, value);
}

bool OpenHurricane::nozzleRegion::distributing(realGeometryArray<cellMesh> &cellQ,
                                               std::string &value) const {
    if (!isOptionSet()) {
        PLWarning("Option unset for mark region: %d", id());
        return false;
    }

    if (!pressPatchedFlag_) {
        auto &mesh = cellQ.mesh();
        const auto varIter = mesh.table().find("p");
        if (varIter != mesh.table().end()) {
            auto ob = varIter->second;
            auto &p = static_cast<realGeometryArray<cellMesh> &>(*ob);
            calPressure(p);
        }
    }

    return true;
}

bool OpenHurricane::nozzleRegion::distributing(vectorGeometryArray<cellMesh> &cellQ,
                                               std::string &value) const {
    if (!isOptionSet()) {
        PLWarning("Option unset for mark region: %d", id());
        return false;
    }

    if (!pressPatchedFlag_) {
        auto &mesh = cellQ.mesh();
        const auto varIter = mesh.table().find("p");
        if (varIter != mesh.table().end()) {
            auto ob = varIter->second;
            auto &p = static_cast<realGeometryArray<cellMesh> &>(*ob);
            calPressure(p);
        }
    }
    return true;
}
