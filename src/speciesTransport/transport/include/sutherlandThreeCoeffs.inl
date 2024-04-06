#include "sutherlandThreeCoeffs.hpp"
/*!
 * \file sutherlandThreeCoeffs.inl
 * \brief The In-Line functions of transport properties by sutherland three coeffs.
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

inline OpenHurricane::sutherlandThreeCoeffs::sutherlandThreeCoeffs(const speciesList &sp,
                                                               const integer index, const real mu0,
                                                               const real T0, const real S,
                                                               const real Prl)
    : transport(sp, index, Prl), mu0_(mu0), T0_(T0), S_(S) {}

inline OpenHurricane::sutherlandThreeCoeffs::sutherlandThreeCoeffs(const speciesList &sp,
                                                               const integer index,
                                                               const controller &cont)
    : transport(sp, index, cont), mu0_(cont.findOrDefault<real>("mu0", 1.716e-5)),
      T0_(cont.findOrDefault<real>("T0", 273.11)), S_(cont.findOrDefault<real>("S", 110.56)) {}

inline OpenHurricane::sutherlandThreeCoeffs::sutherlandThreeCoeffs(const sutherlandThreeCoeffs &tra)
    : transport(tra), mu0_(tra.mu0_), T0_(tra.T0_), S_(tra.S_) {}

inline OpenHurricane::sutherlandThreeCoeffs::sutherlandThreeCoeffs(const sutherlandThreeCoeffs &tra,
                                                               const speciesList &sp)
    : transport(tra, sp), mu0_(tra.mu0_), T0_(tra.T0_), S_(tra.S_) {}

inline OpenHurricane::sutherlandThreeCoeffs::~sutherlandThreeCoeffs() noexcept {}

hur_nodiscard inline OpenHurricane::real OpenHurricane::sutherlandThreeCoeffs::mu(const real p,
                                                                          const real T) const {
    real tt = T / T0_;
    return mu0_ * tt * sqrt(tt) * (T0_ + S_) / (T + S_);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::sutherlandThreeCoeffs::kappa(const real p, const real T, const real cpi) const {
    return mu(p, T) * cpi / Pr();
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::sutherlandThreeCoeffs::kappa(const real p,
                                                                             const real T,
                                                                             const real mui,
                                                                             const real cpi) const {
    return mui * cpi / Pr();
}
