#include "constMuKappa.hpp"
/*!
 * \file constMuKappa.inl
 * \brief The In-Line functions of transport properties by constMuKappa.
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

inline OpenHurricane::constMuKappa::constMuKappa(const speciesList &sp, const integer index,
                                             const real muc, const real kappac, const real Prl)
    : transport(sp, index, Prl), muC_(muc), kappaC_(kappac) {}

inline OpenHurricane::constMuKappa::constMuKappa(const speciesList &sp, const integer index,
                                             const controller &cont)
    : transport(sp, index, cont), muC_(cont.findOrDefault<real>("mu", 1.983e-5)),
      kappaC_(cont.findOrDefault<real>("kappa", 0.0242)) {}

inline OpenHurricane::constMuKappa::constMuKappa(const constMuKappa &tra)
    : transport(tra), muC_(tra.muC_), kappaC_(tra.kappaC_) {}

inline OpenHurricane::constMuKappa::constMuKappa(const constMuKappa &tra, const speciesList &sp)
    : transport(tra, sp), muC_(tra.muC_), kappaC_(tra.kappaC_) {}

inline OpenHurricane::constMuKappa::~constMuKappa() noexcept {}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constMuKappa::mu(const real p, const real T) const {
    return muC_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::constMuKappa::kappa(const real p, const real T,
                                                                    const real cpi) const {
    return kappaC_;
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::constMuKappa::kappa(const real p, const real T, const real mui, const real cpi) const {
    return kappaC_;
}
