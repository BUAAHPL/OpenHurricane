/*!
 * \file ArrheniusThirdBodyRR.cpp
 * \brief Main subroutines for third-body Arrhenius reaction rate.
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

#include "ArrheniusThirdBodyRR.hpp"
namespace OpenHurricane {
    createClassName(ArrheniusThirdBodyRR);
         registerObjFty(reactionRateTypes, ArrheniusThirdBodyRR, controller);
} // namespace OpenHurricane

OpenHurricane::ArrheniusThirdBodyRR::ArrheniusThirdBodyRR(const speciesList &sp, const controller &cont)
    : ArrheniusRR(sp, cont), thirdBodyEff_(sp, cont), M_(0.0) {}

hur_nodiscard OpenHurricane::real OpenHurricane::ArrheniusThirdBodyRR::k(const real p, const real T,
                                                                 const realArray &c) const {
    M_ = thirdBodyEff_.M(c);
    return (M_ * ArrheniusRR::k(p, T, c));
}

void OpenHurricane::ArrheniusThirdBodyRR::gamDGamDCi(const real P, const real T, const realArray &c,
                                                 realArray &gdgdci) const {
#ifdef HUR_DEBUG
    if (c.size() != thirdBodyEff_.size()) {
        LFatal("Size not equal");
    }
#endif // HUR_DEBUG
    if (gdgdci.size() != c.size()) {
        gdgdci.resize(c.size());
    }
    const auto MM = thirdBodyEff_.M(c);
    for (integer i = 0; i < thirdBodyEff_.size(); ++i) {
        gdgdci[i] = thirdBodyEff_[i] / max(MM, tiny);
    }
}
