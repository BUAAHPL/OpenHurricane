/*!
 * \file thirdUnsMUSCL.cpp
 * \brief The subroutines and functions of third order unstructured MUSCL interpolation.
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

#include "thirdUnsMUSCL.hpp"

namespace OpenHurricane {
    createClassNameStr(thirdUnsMUSCL, "thirdUnsMUSCL");
    registerObjFty(reconstruction, thirdUnsMUSCL, controller);
} // namespace OpenHurricane

OpenHurricane::vector OpenHurricane::thirdUnsMUSCL::beta(const vector &r) const {
    vector b;
    for (int i = 0; i < vector::nElements_; ++i) {
        b[i] = beta(r[i]);
    }
    return b;
}

OpenHurricane::vector OpenHurricane::thirdUnsMUSCL::phi(const vector &r) const {
    vector b;
    for (int i = 0; i < vector::nElements_; ++i) {
        b[i] = phi(r[i]);
    }
    return b;
}

OpenHurricane::thirdUnsMUSCL::thirdUnsMUSCL() : reconstruction() {}

OpenHurricane::thirdUnsMUSCL::thirdUnsMUSCL(const controller &cont) : reconstruction(cont) {}

OpenHurricane::thirdUnsMUSCL::~thirdUnsMUSCL() noexcept {}

void OpenHurricane::thirdUnsMUSCL::UnsMUSCLRec(const real QCL, const real QCR, const vector &gradCL,
                                               const vector &gradCR, const vector &dLR, real &UL,
                                               real &UR) const {
    //UI+1 - UI
    real d2 = QCR - QCL;

    //UI - UI-1
    real d1 = 2.0 * (gradCL * dLR) - d2;

    //UI+2 - UI+1
    real d3 = 2.0 * (gradCR * dLR) - d2;

    real rR = (d2 + tiny) / (d3 + tiny);

    real rL = (d2 + tiny) / (d1 + tiny);

    UL = QCL + 0.5 * phi(rL) * d1;
    UR = QCR - 0.5 * phi(rR) * d3;
}

void OpenHurricane::thirdUnsMUSCL::UnsMUSCLRec(
    const vector QCL, const vector QCR, const typename outerProduct<vector, vector>::type &gradCL,
    const typename outerProduct<vector, vector>::type &gradCR, const vector &dLR, vector &UL,
    vector &UR) const {
    //UI+1 - UI
    vector d2 = QCR - QCL;

    //UI - UI-1
    vector d1 = real(2.0) * (gradCL * dLR) - d2;

    //UI+2 - UI+1
    vector d3 = real(2.0) * (gradCR * dLR) - d2;

    vector rR = componentDivide(componentAdd(d2, tiny), componentAdd(d3, tiny));
    vector rL = componentDivide(componentAdd(d2, tiny), componentAdd(d1, tiny));

    UL = QCL + real(0.5) * componentMultiply(phi(rL), d1);
    UR = QCR - real(0.5) * componentMultiply(phi(rR), d3);
}

void OpenHurricane::thirdUnsMUSCL::calcReconstruction(const geometryArray<real, cellMesh> &cellQ,
                                                      const integer faceI, real &ql,
                                                      real &qr) const {
    const runtimeMesh &mesh = cellQ.mesh();
    const faceList &fl = mesh.faces();
    //const auto& cC = mesh.cellCentre();
    const auto &adjoinCCtr = mesh.adjoinCellCtr();

    const integer cl = fl[faceI].leftCell();
    const integer cr = fl[faceI].rightCell();

    //const vector dLR = cC[cr] - cC[cl];
    const vector &dLR = adjoinCCtr[faceI];

    UnsMUSCLRec(cellQ[cl], cellQ[cr], cellQ.grad()[cl], cellQ.grad()[cr], dLR, ql, qr);
}

void OpenHurricane::thirdUnsMUSCL::calcReconstruction(const geometryArray<vector, cellMesh> &cellQ,
                                                      const integer faceI, vector &ql,
                                                      vector &qr) const {
    const runtimeMesh &mesh = cellQ.mesh();
    const faceList &fl = mesh.faces();
    //const auto& cC = mesh.cellCentre();
    const auto &adjoinCCtr = mesh.adjoinCellCtr();

    const integer cl = fl[faceI].leftCell();
    const integer cr = fl[faceI].rightCell();

    //const vector dLR = cC[cr] - cC[cl];
    const vector &dLR = adjoinCCtr[faceI];

    UnsMUSCLRec(cellQ[cl], cellQ[cr], cellQ.grad()[cl], cellQ.grad()[cr], dLR, ql, qr);
}

void OpenHurricane::thirdUnsMUSCL::calcReconstruction(const geometryArray<real, cellMesh> &cellQ,
                                                      const integer faceZoneI, realArray &ql,
                                                      realArray &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();

    const auto &fZ = mesh.faceZones()[faceZoneI];
    //const auto& fC = mesh.faceCentre();
    //const auto& cC = mesh.cellCentre();
    const auto &adjoinCCtr = mesh.adjoinCellCtr();
    const auto &grd = cellQ.grad();
    integer count = 0;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        const integer &cl = fl[fi].leftCell();
        const integer &cr = fl[fi].rightCell();

        //const vector dLR = cC[cr] - cC[cl];
        const vector &dLR = adjoinCCtr[fi];

        //UI+1 - UI
        const real d2 = cellQ[cr] - cellQ[cl];

        //UI - UI-1
        const real d1 = 2.0 * (grd[cl] * dLR) - d2;

        //UI+2 - UI+1
        const real d3 = 2.0 * (grd[cr] * dLR) - d2;

        const real rR = (d2 + tiny) / (d3 + tiny);

        const real rL = (d2 + tiny) / (d1 + tiny);

        ql[count] = cellQ[cl] + 0.5 * phi(rL) * d1;
        qr[count] = cellQ[cr] - 0.5 * phi(rR) * d3;
        count++;
    }
}

void OpenHurricane::thirdUnsMUSCL::calcReconstruction(const geometryArray<vector, cellMesh> &cellQ,
                                                      const integer faceZoneI, vectorArray &ql,
                                                      vectorArray &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();

    const auto &fZ = mesh.faceZones()[faceZoneI];
    const auto &fC = mesh.faceCentre();
    //const auto& cC = mesh.cellCentre();
    const auto &adjoinCCtr = mesh.adjoinCellCtr();
    const auto &grd = cellQ.grad();
    integer count = 0;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        const integer &cl = fl[fi].leftCell();
        const integer &cr = fl[fi].rightCell();

        //const vector dLR = cC[cr] - cC[cl];
        const vector &dLR = adjoinCCtr[fi];

        //UI+1 - UI
        const vector d2 = cellQ[cr] - cellQ[cl];

        //UI - UI-1
        const vector d1 = real(2.0) * (grd[cl] * dLR) - d2;

        //UI+2 - UI+1
        const vector d3 = real(2.0) * (grd[cr] * dLR) - d2;

        const vector rR = componentDivide(componentAdd(d2, tiny), componentAdd(d3, tiny));
        const vector rL = componentDivide(componentAdd(d2, tiny), componentAdd(d1, tiny));

        ql[count] = cellQ[cl] + real(0.5) * componentMultiply(phi(rL), d1);
        qr[count] = cellQ[cr] - real(0.5) * componentMultiply(phi(rR), d3);
        count++;
    }
}
