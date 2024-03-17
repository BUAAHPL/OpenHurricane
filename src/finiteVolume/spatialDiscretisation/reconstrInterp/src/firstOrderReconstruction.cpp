/*!
 * \file firstOrderReconstruction.cpp
 * \brief Main subroutines for first order reconstruction.
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

#include "firstOrderReconstruction.hpp"

namespace OpenHurricane {
    createClassNameStr(firstOrderReconstruction, "firstOrder");
    registerObjFty(reconstruction, firstOrderReconstruction, controller);
} // namespace OpenHurricane

OpenHurricane::firstOrderReconstruction::firstOrderReconstruction() : reconstruction() {}

OpenHurricane::firstOrderReconstruction::firstOrderReconstruction(const controller &cont)
    : reconstruction(cont) {}

OpenHurricane::firstOrderReconstruction::~firstOrderReconstruction() noexcept {}

void OpenHurricane::firstOrderReconstruction::calcReconstruction(
    const geometryArray<real, cellMesh> &cellQ, const integer faceI, real &ql, real &qr) const {
    const runtimeMesh &mesh = cellQ.mesh();
    const faceList &fl = mesh.faces();

    const integer &cl = fl[faceI].leftCell();
    const integer &cr = fl[faceI].rightCell();

    firstOrderRec(cellQ[cl], cellQ[cr], ql, qr);
}

void OpenHurricane::firstOrderReconstruction::calcReconstruction(
    const geometryArray<vector, cellMesh> &cellQ, const integer faceI, vector &ql,
    vector &qr) const {
    const runtimeMesh &mesh = cellQ.mesh();
    const faceList &fl = mesh.faces();

    const integer &cl = fl[faceI].leftCell();
    const integer &cr = fl[faceI].rightCell();

    firstOrderRec(cellQ[cl], cellQ[cr], ql, qr);
}

void OpenHurricane::firstOrderReconstruction::calcReconstruction(
    const geometryArray<real, cellMesh> &cellQ, const integer faceZoneI, realArray &ql,
    realArray &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();

    const auto &fZ = mesh.faceZones()[faceZoneI];

    integer count = 0;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        const integer &cl = fl[fi].leftCell();
        const integer &cr = fl[fi].rightCell();
        ql[count] = cellQ[cl];
        qr[count] = cellQ[cr];
        count++;
    }
}

void OpenHurricane::firstOrderReconstruction::calcReconstruction(
    const geometryArray<vector, cellMesh> &cellQ, const integer faceZoneI, vectorArray &ql,
    vectorArray &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();

    const auto &fZ = mesh.faceZones()[faceZoneI];

    integer count = 0;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        const integer &cl = fl[fi].leftCell();
        const integer &cr = fl[fi].rightCell();
        ql[count] = cellQ[cl];
        qr[count] = cellQ[cr];
        count++;
    }
}
