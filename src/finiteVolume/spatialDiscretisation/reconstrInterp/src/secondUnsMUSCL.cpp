/*!
 * \file secondUnsMUSCL.cpp
 * \brief The subroutines and functions of second order unstructured MUSCL interpolation.
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

#include "secondUnsMUSCL.hpp"
#include "vanLeerLimiter.hpp"

namespace OpenHurricane {
    createClassNameStr(secondUnsMUSCL, "secondUnsMUSCL");
    registerObjFty(reconstruction, secondUnsMUSCL, controller);
} // namespace OpenHurricane

OpenHurricane::secondUnsMUSCL::secondUnsMUSCL()
    : reconstruction(), k_(1.0 / 3.0), limiter_(new vanLeerLimiter()) {}

OpenHurricane::secondUnsMUSCL::secondUnsMUSCL(const controller &cont)
    : reconstruction(cont), k_(cont.findOrDefault<real>("k", 1.0 / 3.0)), limiter_(nullptr) {
    limiter_ = limitersForMUSCL::creator(cont);
}

OpenHurricane::secondUnsMUSCL::~secondUnsMUSCL() noexcept {
    limiter_.clear();
}

void OpenHurricane::secondUnsMUSCL::UnsMUSCLRec(const real QCL, const real QCR,
                                                const vector &gradCL, const vector &gradCR,
                                                const vector &dLR, real &UL, real &UR) const {
    //UI+1 - UI
    const real d2 = QCR - QCL;

    //UI - UI-1
    const real d1 = 2.0 * (gradCL * dLR) - d2;

    //UI+2 - UI+1
    const real d3 = 2.0 * (gradCR * dLR) - d2;

    const real rR = (d2 + tiny) / (d3 + tiny);

    const real rL = (d2 + tiny) / (d1 + tiny);

    UR = QCR -
         real(0.25) * (real((1.0 + k_)) * rR * phi(1.0 / rR) + (real(1.0) - k_) * phi(rR)) * d3;
    UL = QCL +
         real(0.25) * (real((1.0 + k_)) * rL * phi(1.0 / rL) + (real(1.0) - k_) * phi(rL)) * d1;
}

void OpenHurricane::secondUnsMUSCL::UnsMUSCLRec(
    const vector QCL, const vector QCR, const typename outerProduct<vector, vector>::type &gradCL,
    const typename outerProduct<vector, vector>::type &gradCR, const vector &dLR, vector &UL,
    vector &UR) const {
    //UI+1 - UI
    const vector d2 = QCR - QCL;

    //UI - UI-1
    const vector d1 = real(2.0) * (gradCL * dLR) - d2;

    //UI+2 - UI+1
    const vector d3 = real(2.0) * (gradCR * dLR) - d2;

    const vector rR = componentDivide(componentAdd(d2, tiny), componentAdd(d3, tiny));
    const vector rL = componentDivide(componentAdd(d2, tiny), componentAdd(d1, tiny));

    UR =
        QCR - real(0.25) * componentMultiply(((real(1) + k_) * componentMultiply(rR, phi(inv(rR))) +
                                              (real(1.0) - k_) * phi(rR)),
                                             d3);

    UL =
        QCL + real(0.25) * componentMultiply(((real(1) + k_) * componentMultiply(rL, phi(inv(rL))) +
                                              (real(1.0) - k_) * phi(rL)),
                                             d1);
}

void OpenHurricane::secondUnsMUSCL::calcReconstruction(const geometryArray<real, cellMesh> &cellQ,
                                                       const integer faceI, real &ql,
                                                       real &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();
    //const auto& cC = mesh.cellCentre();
    const auto &adjoinCCtr = mesh.adjoinCellCtr();

    const auto &cl = fl[faceI].leftCell();
    const auto &cr = fl[faceI].rightCell();

    //const vector dLR = cC[cr] - cC[cl];
    const vector &dLR = adjoinCCtr[faceI];

    UnsMUSCLRec(cellQ[cl], cellQ[cr], cellQ.grad()[cl], cellQ.grad()[cr], dLR, ql, qr);
}

void OpenHurricane::secondUnsMUSCL::calcReconstruction(const geometryArray<vector, cellMesh> &cellQ,
                                                       const integer faceI, vector &ql,
                                                       vector &qr) const {
    const auto &mesh = cellQ.mesh();
    const auto &fl = mesh.faces();
    //const auto& cC = mesh.cellCentre();
    const auto &adjoinCCtr = mesh.adjoinCellCtr();

    const auto &cl = fl[faceI].leftCell();
    const auto &cr = fl[faceI].rightCell();

    //const vector dLR = cC[cr] - cC[cl];
    const vector &dLR = adjoinCCtr[faceI];

    UnsMUSCLRec(cellQ[cl], cellQ[cr], cellQ.grad()[cl], cellQ.grad()[cr], dLR, ql, qr);
}

void OpenHurricane::secondUnsMUSCL::calcReconstruction(const geometryArray<real, cellMesh> &cellQ,
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
        const real d3 = 2.0 * (grd[cl] * dLR) - d2;

        const real rR = (d2 + tiny) / (d3 + tiny);

        const real rL = (d2 + tiny) / (d1 + tiny);

        qr[count] =
            cellQ[cr] -
            real(0.25) * (real((1.0 + k_)) * rR * phi(1.0 / rR) + (real(1.0) - k_) * phi(rR)) * d3;
        ql[count] =
            cellQ[cl] +
            real(0.25) * (real((1.0 + k_)) * rL * phi(1.0 / rL) + (real(1.0) - k_) * phi(rL)) * d1;

        count++;
    }
}

void OpenHurricane::secondUnsMUSCL::calcReconstruction(const geometryArray<vector, cellMesh> &cellQ,
                                                       const integer faceZoneI, vectorArray &ql,
                                                       vectorArray &qr) const {
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
        const vector d2 = cellQ[cr] - cellQ[cl];

        //UI - UI-1
        const vector d1 = real(2.0) * (grd[cl] * dLR) - d2;

        //UI+2 - UI+1
        const vector d3 = real(2.0) * (grd[cl] * dLR) - d2;

        const vector rR = componentDivide(componentAdd(d2, tiny), componentAdd(d3, tiny));
        const vector rL = componentDivide(componentAdd(d2, tiny), componentAdd(d1, tiny));

        qr[count] =
            cellQ[cr] -
            real(0.25) * componentMultiply(((real(1) + k_) * componentMultiply(rR, phi(inv(rR))) +
                                            (real(1.0) - k_) * phi(rR)),
                                           d3);

        ql[count] =
            cellQ[cl] +
            real(0.25) * componentMultiply(((real(1) + k_) * componentMultiply(rL, phi(inv(rL))) +
                                            (real(1.0) - k_) * phi(rL)),
                                           d1);

        count++;
    }
}
