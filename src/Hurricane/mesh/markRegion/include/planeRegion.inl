#include "planeRegion.hpp"
/*!
 * \file planeRegion.inl
 * \brief In-Line subroutines of the <i>planeRegion.hpp</i> file.
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

inline OpenHurricane::vector OpenHurricane::planeRegion::normal() const noexcept {
    return normal_;
}

inline OpenHurricane::vector OpenHurricane::planeRegion::planePoint() const noexcept {
    return planePoint_;
}

template <class Type, class meshType>
inline bool OpenHurricane::planeRegion::patch(geometryArray<Type, meshType> &cellQ, Type &value) const {
    const auto &mesh = cellQ.mesh();

    if (!isOptionSet()) {
        PLWarning("Option unset for mark region: %d", id());
        return false;
    }

    vector sc = normal_.normalized();
    vector pp = planePoint_;

    const auto &cC = mesh.cellCentre();

    for (integer n = 0; n < mesh.nCells(); ++n) {
        real proj = (cC[n] - pp) * sc;

        if (isInside()) {
            if (proj > 0.0)
                continue;
        } else {
            if (proj <= 0.0)
                continue;
        }
        cellQ[n] = value;
    }

    return true;
}
