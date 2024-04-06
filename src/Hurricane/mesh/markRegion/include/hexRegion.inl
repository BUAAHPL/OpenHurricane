#include "hexRegion.hpp"
/*!
 * \file hexRegion.inl
 * \brief In-Line subroutines of the <i>hexRegion.hpp</i> file.
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

hur_nodiscard inline OpenHurricane::real OpenHurricane::hexRegion::xmax() const noexcept {
    return xmax_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::hexRegion::xmin() const noexcept {
    return xmin_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::hexRegion::ymax() const noexcept {
    return ymax_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::hexRegion::ymin() const noexcept {
    return ymin_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::hexRegion::zmax() const noexcept {
    return zmax_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::hexRegion::zmin() const noexcept {
    return zmin_;
}

template <class Type, class meshType>
inline bool OpenHurricane::hexRegion::patch(geometryArray<Type, meshType> &cellQ, Type &value) const {
    const auto &mesh = cellQ.mesh();

    if (!isOptionSet()) {
        PLWarning("Option unset for mark region: %d",id());
        return false;
    }

    real xma = xmax_;
    real xmi = xmin_;
    real yma = ymax_;
    real ymi = ymin_;
    real zma = zmax_;
    real zmi = zmin_;

    const auto &cC = mesh.cellCentre();

    for (integer n = 0; n < mesh.nCells(); ++n) {
        real xx = cC[n].x();
        real yy = cC[n].y();
        real zz = cC[n].z();

        if (isInside()) {
            if (xx < xmi || xx > xma)
                continue;
            if (yy < ymi || yy > yma)
                continue;
            if (zz < zmi || zz > zma)
                continue;
        } else {
            if (xx >= xmi && xx <= xma && yy >= ymi && yy <= yma && zz >= zmi && zz <= zma)
                continue;
        }
        cellQ[n] = value;
    }

    return true;
}
