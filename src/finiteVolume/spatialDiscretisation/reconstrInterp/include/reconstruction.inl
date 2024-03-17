﻿/*!
 * \file reconstruction.inl
 * \brief In-Line subroutines of the <i>reconstruction.hpp</i> file.
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
#pragma once

inline OpenHurricane::gradient &OpenHurricane::reconstruction::grad() const {
    if (!gradPtr_) {
        LFatal("Attempt to access a null pointer for gradient method class");
    }
    return *gradPtr_;
}

inline void OpenHurricane::reconstruction::calcGrad(geometryArray<real, cellMesh> &cellQ) const {
    gradPtr_->grad(cellQ);
}

inline void OpenHurricane::reconstruction::calcGrad(geometryArray<vector, cellMesh> &cellQ) const {
    gradPtr_->grad(cellQ);
}

inline void
OpenHurricane::reconstruction::calcGrad(PtrList<geometryArray<real, cellMesh>> &cellQList) const {
    for (integer i = 0; i < cellQList.size(); ++i) {
        calcGrad(cellQList[i]);
    }
}

inline void
OpenHurricane::reconstruction::calcGrad(PtrList<geometryArray<vector, cellMesh>> &cellQList) const {
    for (integer i = 0; i < cellQList.size(); ++i) {
        calcGrad(cellQList[i]);
    }
}

inline void OpenHurricane::reconstruction::calcLimiter(geometryArray<real, cellMesh> &cellQ) const {
    // Nothing to do
}

inline void
OpenHurricane::reconstruction::calcLimiter(geometryArray<vector, cellMesh> &cellQ) const {
    // Nothing to do
}

template <class Type>
inline void
OpenHurricane::reconstruction::calcLimiter(PtrList<geometryArray<Type, cellMesh>> &cellQ) const {
    for (integer i = 0; i < cellQ.size(); ++i) {
        calcLimiter(cellQ[i]);
    }
}
