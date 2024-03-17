/*!
 * \file limitersForLinear.inl
 * \brief In-Line subroutines of the <i>limitersForLinear.hpp</i> file.
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

template <class Type>
inline void OpenHurricane::limitersForLinear::limiters(
    const geometryArray<Type, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &gradQ,
    geometryArray<Type, cellMesh> &_limiters) const {
    calcLimiters(cellQ, gradQ, _limiters);
    updateBoundaryField(_limiters);
}

template <class Type>
inline void OpenHurricane::limitersForLinear::limiters(geometryArray<Type, cellMesh> &cellQ) const {
    calcLimiters(cellQ, cellQ.grad(), cellQ.limiter());
    updateBoundaryField(cellQ.limiter());
}

template <class Type>
inline OpenHurricane::geometryArray<Type, OpenHurricane::cellMesh>
OpenHurricane::limitersForLinear::limiters(
    const geometryArray<Type, cellMesh> &cellQ,
    const geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &gradQ) const {
    geometryArray<Type, cellMesh> limi(cellQ.size());

    limiters(cellQ, gradQ, limi);

    return limi;
}
