/*!
 * \file JanevRR.hpp
 * \brief Main subroutines for Janev reaction rate.
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

#include "JanevRR.hpp"

namespace OpenHurricane {
    createClassName(JanevRR);
         registerObjFty(reactionRateTypes, JanevRR, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::JanevRR::k(const real p, const real T,
                                                    const realArray &c) const {
    real lta = A_;

    if (mag(beta_) > veryTiny) {
        lta *= std::pow(T, beta_);
    }
    real expArg = 0.0;

    if (mag(Ta_) > veryTiny) {
        expArg -= Ta_ / T;
    }

    real lnT = std::log(T);

    for (int n = 0; n < nb_; n++) {
        expArg += b_[n] * std::pow(lnT, n);
    }

    lta *= std::exp(expArg);
    return lta;
}

hur_nodiscard OpenHurricane::real OpenHurricane::JanevRR::DkDT(const real kj, const real p, const real T,
                                            const realArray &c) const {
    real lta = A_;

    if (mag(beta_) > veryTiny) {
        lta *= pow(T, beta_);
    }

    real expArg = 0;

    if (mag(Ta_) > veryTiny) {
        expArg -= Ta_ / T;
    }

    real lnT = log(T);

    for (int n = 0; n < nb_; n++) {
        expArg += b_[n] * pow(lnT, n);
    }

    real deriv = b_[1];

    for (int n = 2; n < nb_; n++) {
        deriv += n * b_[n] * pow(lnT, n - 1);
    }

    lta *= exp(expArg);

    return lta * (beta_ + Ta_ / T + deriv) / T;
}