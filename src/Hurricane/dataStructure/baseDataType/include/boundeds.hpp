/*!
 * \file boundeds.hpp
 * \brief Headers of boundeds data.
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

#include "bounded.hpp"
#include "vectors.hpp"

namespace OpenHurricane {
    /**
     * \brief Bounded value.
     */
    using realBounded = bounded<real>;

    /**
     * \brief Bounded value.
     */
    using vectorBounded = bounded<vector>;

    template <> void realBounded::checkBoundValue() const;

    template <> void vectorBounded::checkBoundValue() const;

    template <>
    hur_nodiscard real OpenHurricane::bounded<real>::bounding(const real &v, const real &limit) const;

    template <>
    hur_nodiscard vector OpenHurricane::bounded<vector>::bounding(const vector &v,
                                                              const vector &limit) const;
} // namespace OpenHurricane
