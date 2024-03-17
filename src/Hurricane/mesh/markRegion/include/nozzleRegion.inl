#include "nozzleRegion.hpp"
/*!
 * \file nozzleRegion.inl
 * \brief In-Line subroutines of the <i>nozzleRegion.hpp</i> file.
 * \author Peng Jian
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

inline const OpenHurricane::vector &OpenHurricane::nozzleRegion::Aexit() const noexcept {
    return Aexit_;
}

hur_nodiscard inline const OpenHurricane::vector &OpenHurricane::nozzleRegion::Ahead() const noexcept {
    return Ahead_;
}

hur_nodiscard inline const OpenHurricane::vector &OpenHurricane::nozzleRegion::Athroat() const noexcept {
    return Athroat_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::nozzleRegion::Dhead() const noexcept {
    return Dhead_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::nozzleRegion::Dexit() const noexcept {
    return Dexit_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::nozzleRegion::Dthroat() const noexcept {
    return Dthroat_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::nozzleRegion::gammaConst() const noexcept {
    return gammaConst_;
}

hur_nodiscard inline const std::string &OpenHurricane::nozzleRegion::bcfName() const noexcept {
    return boundaryFaceName_;
}

template <class Type, class meshType>
inline bool OpenHurricane::nozzleRegion::patch(geometryArray<Type, meshType> &cellQ,
                                           Type &value) const {
    const auto &mesh = cellQ.mesh();

    if (!isOptionSet()) {
        PLWarning("Option unset for mark region: %d", id());
        return false;
    }

    vector Ae = Aexit_;
    vector Ah = Ahead_;
    vector At = Athroat_;
    real De = Dexit_;
    real Dh = Dhead_;
    real Dt = Dthroat_;

    real AhAe = dist(Ah, Ae);

    const auto &cC = mesh.cellCentre();

    for (integer n = 0; n < mesh.nCells(); ++n) {
        const auto &C = cC[n];

        vector d = ((Ae - Ah) ^ (C - Ah)) / AhAe;
        real nLen = d.magnitude();

        // AC * AB
        real CABA = (C - Ah) * (Ae - Ah);

        // BC * BA
        real CBAB = (C - Ae) * (Ah - Ae);

        if (isInside()) {
            if (nLen > De)
                continue;
            if (CABA < 0.0 || CBAB < 0.0)
                continue;
        } else {
            if (nLen <= De && CABA > 0.0 && CBAB >= 0.0)
                continue;
        }

        cellQ[n] = value;
    }

    return true;
}