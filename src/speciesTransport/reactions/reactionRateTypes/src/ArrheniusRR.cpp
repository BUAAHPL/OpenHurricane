/*!
 * \file ArrheniusRR.cpp
 * \brief Main subroutines for Arrhenius reaction rate.
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

#include "ArrheniusRR.hpp"

namespace OpenHurricane {
    createClassName(ArrheniusRR);
         registerObjFty(reactionRateTypes, ArrheniusRR, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::ArrheniusRR::k(const real p, const real T,
                                                      const realArray &c) const {
    real k_ = A_;

    if (mag(beta_) > veryTiny) {
        k_ *= std::pow(T, beta_);
    }
    if (mag(Ta_) > veryTiny) {
        k_ *= exp(-Ta_ / T);
    }
    return k_;
}