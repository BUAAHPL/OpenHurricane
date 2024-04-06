/*!
 * \file powerSeriesRR.hpp
 * \brief Main subroutines for power series reaction rate.
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

#include "powerSeriesRR.hpp"

namespace OpenHurricane {
    createClassName(powerSeriesRR);
         registerObjFty(reactionRateTypes, powerSeriesRR, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::powerSeriesRR::k(const real p, const real T,
                                                          const realArray &c) const {
    real lta = A_;

    if (mag(beta_) > veryTiny) {
        lta *= pow(T, beta_);
    }
    real expArg = 0.0;

    for (int n = 0; n < nb_; n++) {
        expArg += b_[n] / pow(T, n + 1);
    }

    lta *= exp(expArg);
    return lta;
}

hur_nodiscard OpenHurricane::real OpenHurricane::powerSeriesRR::DkDT(const real kj, const real p,
                                                             const real T,
                                                             const realArray &c) const {
    real lta = A_;

    if (mag(beta_) > veryTiny) {
        lta *= pow(T, beta_);
    }

    real expArg = 0;
    real deriv = 0;

    for (int n = 0; n < nb_; n++) {
        real cT = b_[n] / pow(T, n + 1);
        expArg += cT;
        deriv -= (n + 1) * cT;
    }

    lta *= exp(expArg);

    return lta * (beta_ + deriv) / T;
}