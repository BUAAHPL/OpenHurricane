/*!
 * \file LandauTellerRR.hpp
 * \brief Main subroutines for Landau-Teller reaction rate.
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

#include "LandauTellerRR.hpp"

namespace OpenHurricane {
    createClassName(LandauTellerRR);
         registerObjFty(reactionRateTypes, LandauTellerRR, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::LandauTellerRR::k(const real p, const real T,
                                                           const realArray &c) const {
    real lta = A_;

    if (mag(beta_) > veryTiny) {
        lta *= pow(T, beta_);
    }

    real expAgr = 0.0;
    if (mag(Ta_) > veryTiny) {
        expAgr -= Ta_ / T;
    }

    if (mag(B_) > veryTiny) {
        expAgr += B_ / pow(T, real(1.0 / 3.0));
    }

    if (mag(C_) > veryTiny) {
        expAgr += C_ / pow(T, real(2.0 / 3.0));
    }

    if (mag(expAgr) > veryTiny) {
        lta *= exp(expAgr);
    }

    return lta;
}

hur_nodiscard OpenHurricane::real OpenHurricane::LandauTellerRR::DkDT(const real kj, const real p,
                                                              const real T,
                                                              const realArray &c) const {
    real lta = A_;

    if (mag(beta_) > veryTiny) {
        lta *= pow(T, beta_);
    }

    real expArg = 0;
    real deriv = 0;

    if (mag(Ta_) > veryTiny) {
        real TaT = Ta_ / T;
        expArg -= TaT;
        deriv += TaT;
    }

    if (mag(B_) > veryTiny) {
        real BT = B_ / cbrt(T);
        expArg += BT;
        deriv -= BT / 3;
    }

    if (mag(C_) > veryTiny) {
        real CT = C_ / pow(T, real(2.0 / 3.0));
        expArg += CT;
        deriv -= 2 * CT / 3;
    }

    if (mag(expArg) > veryTiny) {
        lta *= exp(expArg);
    }

    return lta * (beta_ + deriv) / T;
}