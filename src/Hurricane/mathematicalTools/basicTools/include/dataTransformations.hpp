/*!
 * \file dataTransformations.hpp
 * \brief Header of data transformations.
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
#include "OpenHurricane.hpp"

namespace OpenHurricane {
    namespace dataTransformations {
        hur_nodiscard constexpr inline real maxMin(const real x, const real xmax,
                                                   const real xmin) noexcept {
#ifdef HUR_DEBUG
            if (xmax <= xmin) {
                LFatal("The value of xmax: %e must be larger than xmin: %e", xmax, xmin);
            }
#endif // HUR_DEBUG

            return (x - xmin) / (xmax - xmin);
        }

        hur_nodiscard inline real BoxCox(const real x, const real lambda) {
            if (lambda == 0) {
                return log(x);
            } else if (lambda == 1) {
                return x - 1.0;
            } else if (lambda == 2) {
                return (sqr(x) - 1.0) / 2.0;
            } else {
                return (pow(x, lambda) - 1.0) / lambda;
            }
        }

    } // namespace dataTransformations
} //  namespace OpenHurricane
