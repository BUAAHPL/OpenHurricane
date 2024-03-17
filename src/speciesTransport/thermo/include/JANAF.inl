/*!
 * \file JANAF.inl
 * \brief The In-Line functions of thermo properties by JANAF.
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

#include "JANAF.hpp"
#include "species.hpp"

hur_nodiscard inline const typename OpenHurricane::JANAF::coeffArray &
OpenHurricane::JANAF::coeffs(const real T) const noexcept {
    if (T < TCommon_) {
        return lowCpCoeffs_;
    } else {
        return highCpCoeffs_;
    }
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::JANAF::limit(const real T) const noexcept {
    if (T < TLow_ || T > THigh_) {
#ifdef _FULLDEBUG
        Info << " Exceed temperature range of JANAF thermo: " << TLow_ << " ~ " << THigh_
             << ";  T = " << T << std::endl;
#endif // _FULLDEBUG
        return min(max(T, TLow_), THigh_);
    } else {
        return T;
    }
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::JANAF::Tlow() const noexcept {
    return TLow_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::JANAF::Thigh() const noexcept {
    return THigh_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::JANAF::Tcommon() const noexcept {
    return TCommon_;
}

hur_nodiscard inline const typename OpenHurricane::JANAF::coeffArray &
OpenHurricane::JANAF::highCpCoeffs() const noexcept {
    return highCpCoeffs_;
}

hur_nodiscard inline const typename OpenHurricane::JANAF::coeffArray &
OpenHurricane::JANAF::lowCpCoeffs() const noexcept {
    return lowCpCoeffs_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::JANAF::Dha0DT(const real T) const noexcept {
    return cp0(T);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::JANAF::Ds0DT(const real T) const noexcept {
    const coeffArray &a = coeffs(T);
    return (a[0] / max(T, tiny) + a[1] + (a[2] + (a[3] + a[4] * T) * T) * T);
}
