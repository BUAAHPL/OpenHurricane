/*!
 * \file constCv.inl
 * \brief The In-Line functions of thermo properties by constCv.
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

#include "constCv.hpp"
#include "species.hpp"

inline OpenHurricane::constCv::constCv(const equationOfState &st, const real cv, const real hc,
                                   const integer id)
    : thermo(st, id), cv_(cv), hc_(hc) {}

inline OpenHurricane::constCv::constCv(const constCv &jt) : thermo(jt), cv_(jt.cv_), hc_(jt.hc_) {}

inline OpenHurricane::constCv::constCv(const constCv &jt, const string &name)
    : thermo(jt), cv_(jt.cv_), hc_(jt.hc_) {}

hur_nodiscard OpenHurricane::real OpenHurricane::constCv::cp0(const real T) const noexcept {
    return cv0(T) + eos().species()[speciesIndex_].Ri();
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constCv::cv0(const real T) const noexcept {
    return cv_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constCv::ha0(const real T) const noexcept {
    return cp0(T) * T + hc_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constCv::Dha0DT(const real T) const noexcept {
    return cp0(T);
}

hur_nodiscard OpenHurricane::real OpenHurricane::constCv::hc() const noexcept {
    return hc_;
}

hur_nodiscard OpenHurricane::real OpenHurricane::constCv::s0(const real T) const {
    return cp0(T) * log(max(T, tiny));
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constCv::Ds0DT(const real T) const noexcept {
    return cp0(T) / max(T, tiny);
}
