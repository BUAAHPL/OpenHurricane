/*!
 * \file constCp.inl
 * \brief The In-Line functions of thermo properties by constCp.
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

#include "constCp.hpp"
#include "species.hpp"

inline OpenHurricane::constCp::constCp(const equationOfState &st, const real cp, const real hc,
                                   const integer id)
    : thermo(st, id), cp_(cp), hc_(hc) {}

inline OpenHurricane::constCp::constCp(const constCp &jt) : thermo(jt), cp_(jt.cp_), hc_(jt.hc_) {}

inline OpenHurricane::constCp::constCp(const constCp &jt, const string &name)
    : thermo(jt), cp_(jt.cp_), hc_(jt.hc_) {}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constCp::cp0(const real T) const noexcept{
    return cp_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constCp::ha0(const real T) const noexcept {
    return cp0(T) * T + hc_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constCp::Dha0DT(const real T) const noexcept {
    return cp0(T);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constCp::hc() const noexcept {
    return hc_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constCp::s0(const real T) const {
    return cp0(T) * log(max(T, tiny));
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constCp::Ds0DT(const real T) const noexcept {
    return cp0(T) / max(T, tiny);
}