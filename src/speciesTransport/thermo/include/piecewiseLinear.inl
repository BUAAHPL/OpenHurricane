/*!
 * \file piecewiseLinear.inl
 * \brief The In-Line functions of thermo properties by piecewise-linear.
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

#include "piecewiseLinear.hpp"
#include "species.hpp"

inline OpenHurricane::piecewiseLinear::piecewiseLinear(const equationOfState &st, const realArray cpl,
                                                   const realArray Tl, const real hc,
                                                   const integer id)
    : thermo(st, id), cp_(cpl), T_(Tl), hc_(hc) {}

inline OpenHurricane::piecewiseLinear::piecewiseLinear(const piecewiseLinear &jt)
    : thermo(jt), cp_(jt.cp_), T_(jt.T_), hc_(jt.hc_) {}

inline OpenHurricane::piecewiseLinear::piecewiseLinear(const piecewiseLinear &jt, const string &name)
    : thermo(jt), cp_(jt.cp_), T_(jt.T_), hc_(jt.hc_) {}

hur_nodiscard inline OpenHurricane::real OpenHurricane::piecewiseLinear::cp0(const real T) const noexcept{
    if (T <= T_[0]) {
        return cp_[0];
    } else if (T >= T_.last()) {
        return cp_.last();
    }
    for (integer i = 0; i < T_.size() - 1; ++i) {
        if (T >= T_[i] && T < T_[i + 1]) {
            return cp_[i] + (cp_[i + 1] - cp_[i]) * (T - T_[i]) / (T_[i + 1] - T_[i]);
        }
    }
    return cp_[0];
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::piecewiseLinear::inteCp0dT(const real T1, const real T2) const noexcept {
    return real(0);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::piecewiseLinear::ha0(const real T) const noexcept {
    return cp0(T) * T + hc_;
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::piecewiseLinear::Dha0DT(const real T) const noexcept {
    return cp0(T);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::piecewiseLinear::hc() const noexcept {
    return hc_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::piecewiseLinear::s0(const real T) const {
    return cp0(T) * log(max(T, tiny));
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::piecewiseLinear::Ds0DT(const real T) const noexcept {
    return cp0(T) / max(T, tiny);
}