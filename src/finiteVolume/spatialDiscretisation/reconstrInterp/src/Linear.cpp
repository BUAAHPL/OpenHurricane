/*!
 * \file Linear.cpp
 * \brief The subroutines and functions of piecewise linear reconstruction.
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

#include "Linear.hpp"

namespace OpenHurricane {
    createClassNameStr(Linear, "linear");
    registerObjFty(reconstruction, Linear, controller);
} // namespace OpenHurricane

OpenHurricane::Linear::Linear() : reconstruction(), limiterPtr_(nullptr) {}

OpenHurricane::Linear::Linear(const controller &cont) : reconstruction(cont), limiterPtr_(nullptr) {
    limiterPtr_ = limitersForLinear::creator(cont);
}

OpenHurricane::Linear::~Linear() noexcept {}

inline void OpenHurricane::Linear::calcReconstruction(const geometryArray<real, cellMesh> &cellQ,
                                                      const integer faceI, real &ql,
                                                      real &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();
    //const auto& fC = mesh.faceCentre();
    //const auto& cC = mesh.cellCentre();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    const integer &cl = fl[faceI].leftCell();
    const integer &cr = fl[faceI].rightCell();

    /*const vector rL = fC[faceI] - cC[cl];
    const vector rR = fC[faceI] - cC[cr];*/
    const vector &rL = faceLeftCellCtr[faceI];
    const vector &rR = faceRightCellCtr[faceI];

    auto &grd = cellQ.grad();
    auto &lmt = cellQ.limiter();
    linearRecon(cellQ[cl], cellQ[cr], grd[cl], grd[cr], rL, rR, lmt[cl], lmt[cr], ql, qr);
}

void OpenHurricane::Linear::calcReconstruction(const geometryArray<vector, cellMesh> &cellQ,
                                               const integer faceI, vector &ql, vector &qr) const {
    /*const runtimeMesh& mesh = cellQ.mesh();
    const faceList& fl = mesh.faces();
    const vectorArray& fC = mesh.fC();
    const vectorArray& cC = mesh.cC();*/

    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();
    //const auto& fC = mesh.faceCentre();
    //const auto& cC = mesh.cellCentre();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    const integer &cl = fl[faceI].leftCell();
    const integer &cr = fl[faceI].rightCell();

    //const vector rL = fC[faceI] - cC[cl];
    //const vector rR = fC[faceI] - cC[cr];
    const vector &rL = faceLeftCellCtr[faceI];
    const vector &rR = faceRightCellCtr[faceI];

    auto &grd = cellQ.grad();
    auto &lmt = cellQ.limiter();
    linearRecon(cellQ[cl], cellQ[cr], grd[cl], grd[cr], rL, rR, lmt[cl], lmt[cr], ql, qr);
}

void OpenHurricane::Linear::calcReconstruction(const geometryArray<real, cellMesh> &cellQ,
                                               const integer faceZoneI, realArray &ql,
                                               realArray &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();

    const auto &fZ = mesh.faceZones()[faceZoneI];
    //const auto& fC = mesh.faceCentre();
    //const auto& cC = mesh.cellCentre();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    const auto &grd = cellQ.grad();
    const auto &lmt = cellQ.limiter();
    integer count = 0;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        const integer &cl = fl[fi].leftCell();
        const integer &cr = fl[fi].rightCell();
        //const vector rL = fC[fi] - cC[cl];
        //const vector rR = fC[fi] - cC[cr];
        const vector &rL = faceLeftCellCtr[fi];
        const vector &rR = faceRightCellCtr[fi];
        //ql[count] = cellQ[cl] + componentMultiply(lmt[cl], (grd[cl] * rL));
        //qr[count] = cellQ[cr] + componentMultiply(lmt[cr], (grd[cr] * rR));
        ql[count] = cellQ[cl] + lmt[cl] * (grd[cl] * rL);
        qr[count] = cellQ[cr] + lmt[cr] * (grd[cr] * rR);
        count++;
    }
}

void OpenHurricane::Linear::calcReconstruction(const geometryArray<vector, cellMesh> &cellQ,
                                               const integer faceZoneI, vectorArray &ql,
                                               vectorArray &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();

    const auto &fZ = mesh.faceZones()[faceZoneI];
    //const auto& fC = mesh.faceCentre();
    //const auto& cC = mesh.cellCentre();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    const auto &grd = cellQ.grad();
    const auto &lmt = cellQ.limiter();
    integer count = 0;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        const integer &cl = fl[fi].leftCell();
        const integer &cr = fl[fi].rightCell();
        /*const vector rL = fC[fi] - cC[cl];
        const vector rR = fC[fi] - cC[cr];*/
        const vector &rL = faceLeftCellCtr[fi];
        const vector &rR = faceRightCellCtr[fi];
        ql[count] = cellQ[cl] + componentMultiply(lmt[cl], (grd[cl] * rL));
        qr[count] = cellQ[cr] + componentMultiply(lmt[cr], (grd[cr] * rR));

        count++;
    }
}

void OpenHurricane::Linear::calcReconstructionWithoutLimit(
    const geometryArray<real, cellMesh> &cellQ, const integer faceI, real &ql, real &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    const integer &cl = fl[faceI].leftCell();
    const integer &cr = fl[faceI].rightCell();

    const vector &rL = faceLeftCellCtr[faceI];
    const vector &rR = faceRightCellCtr[faceI];

    auto &grd = cellQ.grad();

    ql = cellQ[cl] + grd[cl] * rL;
    qr = cellQ[cr] + grd[cr] * rR;
}

void OpenHurricane::Linear::calcReconstructionWithoutLimit(
    const geometryArray<vector, cellMesh> &cellQ, const integer faceI, vector &ql,
    vector &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();

    const integer &cl = fl[faceI].leftCell();
    const integer &cr = fl[faceI].rightCell();

    const vector &rL = faceLeftCellCtr[faceI];
    const vector &rR = faceRightCellCtr[faceI];

    auto &grd = cellQ.grad();

    ql = cellQ[cl] + grd[cl] * rL;
    qr = cellQ[cr] + grd[cr] * rR;
}

void OpenHurricane::Linear::calcReconstructionWithoutLimit(
    const geometryArray<real, cellMesh> &cellQ, const integer faceZoneI, realArray &ql,
    realArray &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();

    const auto &fZ = mesh.faceZones()[faceZoneI];
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    const auto &grd = cellQ.grad();
    integer count = 0;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        const integer &cl = fl[fi].leftCell();
        const integer &cr = fl[fi].rightCell();
        const vector &rL = faceLeftCellCtr[fi];
        const vector &rR = faceRightCellCtr[fi];
        ql[count] = cellQ[cl] + grd[cl] * rL;
        qr[count] = cellQ[cr] + grd[cr] * rR;
        count++;
    }
}

void OpenHurricane::Linear::calcReconstructionWithoutLimit(
    const geometryArray<vector, cellMesh> &cellQ, const integer faceZoneI, vectorArray &ql,
    vectorArray &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();

    const auto &fZ = mesh.faceZones()[faceZoneI];
    //const auto& fC = mesh.faceCentre();
    //const auto& cC = mesh.cellCentre();
    const auto &faceLeftCellCtr = mesh.faceCtrToLeftCellCtr();
    const auto &faceRightCellCtr = mesh.faceCtrToRightCellCtr();
    const auto &grd = cellQ.grad();
    integer count = 0;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        const integer &cl = fl[fi].leftCell();
        const integer &cr = fl[fi].rightCell();
        const vector &rL = faceLeftCellCtr[fi];
        const vector &rR = faceRightCellCtr[fi];
        ql[count] = cellQ[cl] + grd[cl] * rL;
        qr[count] = cellQ[cr] + grd[cr] * rR;

        count++;
    }
}
