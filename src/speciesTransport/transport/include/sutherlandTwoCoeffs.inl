#include "sutherlandTwoCoeffs.hpp"
/*!
 * \file sutherlandTwoCoeffs.inl
 * \brief The In-Line functions of transport properties by sutherland two coeffs.
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

inline OpenHurricane::sutherlandTwoCoeffs::sutherlandTwoCoeffs(const speciesList &sp,
                                                           const integer index, const real C1,
                                                           const real C2, const real Prl)
    : transport(sp, index, Prl), C1_(C1), C2_(C2) {}

inline OpenHurricane::sutherlandTwoCoeffs::sutherlandTwoCoeffs(const speciesList &sp,
                                                           const integer index,
                                                           const controller &cont)
    : transport(sp, index, cont), C1_(cont.findOrDefault<real>("C1", 1.458e-6)),
      C2_(cont.findOrDefault<real>("C2", 110.4)) {}

inline OpenHurricane::sutherlandTwoCoeffs::sutherlandTwoCoeffs(const sutherlandTwoCoeffs &tra)
    : transport(tra), C1_(tra.C1_), C2_(tra.C2_) {}

inline OpenHurricane::sutherlandTwoCoeffs::sutherlandTwoCoeffs(const sutherlandTwoCoeffs &tra,
                                                           const speciesList &sp)
    : transport(tra, sp), C1_(tra.C1_), C2_(tra.C2_) {}

inline OpenHurricane::sutherlandTwoCoeffs::~sutherlandTwoCoeffs() noexcept {}

hur_nodiscard inline OpenHurricane::real OpenHurricane::sutherlandTwoCoeffs::mu(const real p,
                                                                        const real T) const {
    return C1_ * T * sqrt(T) / (T + C2_);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::sutherlandTwoCoeffs::kappa(const real p, const real T, const real cpi) const {
    return mu(p, T) * cpi / Pr();
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::sutherlandTwoCoeffs::kappa(const real p,
                                                                           const real T,
                                                                           const real mui,
                                                                           const real cpi) const {
    return mui * cpi / Pr();
}
