/*!
 * \file faceGrad.cpp
 * \brief Main subroutines for gradient on face.
 * \author Rao Sihang
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

#include "faceGrad.hpp"

OpenHurricane::faceGrad::faceGrad() {}

typename OpenHurricane::outerProduct<OpenHurricane::vector, OpenHurricane::real>::type
OpenHurricane::faceGrad::grad(
    const geometryArray<real, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &cellGrad,
    const integer faceI) {
    typedef typename outerProduct<vector, real>::type gradType;
    const auto &mesh = cellQ.mesh();
    const auto &faces = mesh.faces();
    const auto &adjoinCCtr = mesh.adjoinCellCtr();

    const auto &fWgt = mesh.faceWeight();

    const auto &cl = faces[faceI].leftCell();
    const auto &cr = faces[faceI].rightCell();

    // Distance vector pointing from cell centre of cr to cl.
    const vector &drl = adjoinCCtr[faceI];

    const real wl = fWgt[faceI];
    const real wr = real(1.0) - wl;

    const gradType averageGrad = wl * cellGrad[cl] + wr * cellGrad[cr];

    const real dqq = (cellQ[cr] - cellQ[cl]) / drl.magnitude();
    const vector drlUnit = drl.normalized();

    const gradType faceGradI = averageGrad - (averageGrad * drlUnit - dqq) * drlUnit;

    return faceGradI;
}

void OpenHurricane::faceGrad::grad(
    const geometryArray<real, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &cellGrad,
    geometryArray<typename outerProduct<vector, real>::type, faceMesh> &faceGrads) {
    const integer nFaces = cellQ.mesh().nFaces();
    typedef typename outerProduct<vector, real>::type gradType;
    const auto &mesh = cellQ.mesh();
    const auto &faces = mesh.faces();
    const auto &adjoinCCtr = mesh.adjoinCellCtr();

    const auto &fWgt = mesh.faceWeight();

    for (integer faceI = 0; faceI < nFaces; ++faceI) {
        const auto &cl = faces[faceI].leftCell();
        const auto &cr = faces[faceI].rightCell();

        // Distance vector pointing from cell centre of cr to cl.
        const vector &drl = adjoinCCtr[faceI];

        const real wl = fWgt[faceI];
        const real wr = real(1.0) - wl;

        const gradType averageGrad = wl * cellGrad[cl] + wr * cellGrad[cr];

        const real dqq = (cellQ[cr] - cellQ[cl]) / drl.magnitude();
        const vector drlUnit = drl.normalized();

        faceGrads[faceI] = averageGrad - (averageGrad * drlUnit - dqq) * drlUnit;
    }
}

typename OpenHurricane::outerProduct<OpenHurricane::vector, OpenHurricane::vector>::type
OpenHurricane::faceGrad::grad(
    const geometryArray<vector, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &cellGrad,
    const integer faceI) {
    typedef typename outerProduct<vector, vector>::type gradType;
    const auto &mesh = cellQ.mesh();
    const auto &faces = mesh.faces();
    const auto &adjoinCCtr = mesh.adjoinCellCtr();

    const auto &fWgt = mesh.faceWgt();

    const auto &cl = faces[faceI].leftCell();
    const auto &cr = faces[faceI].rightCell();

    // Distance vector pointing from cell centre of cr to cl.
    const vector &drl = adjoinCCtr[faceI];

    const real wl = fWgt[faceI];
    const real wr = real(1.0) - wl;

    const gradType averageGrad = wl * cellGrad[cl] + wr * cellGrad[cr];

    const vector dqq = (cellQ[cr] - cellQ[cl]) / drl.magnitude();
    const vector drlUnit = drl.normalized();

    const gradType faceGradI = averageGrad - ((averageGrad * drlUnit - dqq) & drlUnit);

    return faceGradI;
}

void OpenHurricane::faceGrad::grad(
    const geometryArray<vector, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &cellGrad,
    geometryArray<typename outerProduct<vector, vector>::type, faceMesh> &faceGrads) {
    const integer nFaces = cellQ.mesh().nFaces();
    typedef typename outerProduct<vector, vector>::type gradType;
    const auto &mesh = cellQ.mesh();
    const auto &faces = mesh.faces();
    const auto &fC = mesh.faceCentre();
    const auto &adjoinCCtr = mesh.adjoinCellCtr();

    const auto &fWgt = mesh.faceWeight();

    for (integer faceI = 0; faceI < nFaces; ++faceI) {
        const auto &cl = faces[faceI].leftCell();
        const auto &cr = faces[faceI].rightCell();
        // Distance vector pointing from cell centre of cr to cl.
        const vector &drl = adjoinCCtr[faceI];

        const real wl = fWgt[faceI];
        const real wr = real(1.0) - wl;

        const gradType averageGrad = wl * cellGrad[cl] + wr * cellGrad[cr];

        const vector dqq = (cellQ[cr] - cellQ[cl]) / drl.magnitude();
        const vector drlUnit = drl.normalized();

        faceGrads[faceI] = averageGrad - ((averageGrad * drlUnit - dqq) & drlUnit);
    }
}