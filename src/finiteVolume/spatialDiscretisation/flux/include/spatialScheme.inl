#include "spatialScheme.hpp"
/*!
 * \file spatialScheme.inl
 * \brief In-Line subroutines of the <i>spatialScheme.hpp</i> file.
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
#pragma once

hur_nodiscard inline const OpenHurricane::runtimeMesh &
OpenHurricane::spatialScheme::mesh() const noexcept {
    return mesh_;
}

template <class Type>
inline void
OpenHurricane::spatialScheme::invFluxTemplate(geometryArray<Type, cellMesh> &cellQ) const {
    cellQ.rhs() = Zero;
    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);
    const auto &mesh = cellQ.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();
    const auto &fC = mesh.faceCntr();

    const auto &fA = mesh.faceArea();
    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            continue;
        }
        //Other cell face inviscid flux calculations
        else {
            Array<Type> ql(fzl[fzi].size());
            Array<Type> qr(fzl[fzi].size());
            reconstrPtr_->calcReconstruction(cellQ, fzi, ql, qr);
            integer count = 0;
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();

                Type flux;
                if (rhoFlux_[fi] >= 0.0) {
                    flux = rhoFlux_[fi] * qr[count];
                } else {
                    flux = rhoFlux_[fi] * ql[count];
                }
                //Add the numerical flux to the right-hand-term of both cells sharing the face
                cellQ.rhs()[cl] += flux;
                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= flux;
                }
                count++;
            }
        }
    }
}

template <class Type>
void OpenHurricane::spatialScheme::addInvFlux(geometryArray<Type, cellMesh> &cellQ) {
    objectList_.append(cellQ);
    for (int i = 0; i < cellQ.nElements(); ++i) {
        paramMap_[countParam_][0] = objectList_.size() - 1;
        paramMap_[countParam_][1] = i;
        countParam_++;
    }
}

template <class Type>
inline void OpenHurricane::spatialScheme::grad(geometryArray<Type, cellMesh> &cellQ) const {
    reconstrPtr_->calcGrad(cellQ);
}

hur_nodiscard inline const OpenHurricane::flowModel &
OpenHurricane::spatialScheme::flows() const noexcept {
    return flows_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::spatialScheme::MaInf() const noexcept {
    return flows_.refValues().Ma();
}
