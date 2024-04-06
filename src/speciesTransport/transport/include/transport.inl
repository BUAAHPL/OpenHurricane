#include "transport.hpp"
/*!
 * \file transport.inl
 * \brief The In-Line functions of transport properties.
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

hur_nodiscard inline const OpenHurricane::speciesList &OpenHurricane::transport::species() const noexcept {
    return species_;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::transport::index() const noexcept {
    return index_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::transport::Pr() const noexcept {
    return Pr_;
}

inline OpenHurricane::real OpenHurricane::transport::Di(const real p, const real T, const realArray &xi,
                                                const PtrList<transport> &tranPtr) const {
    return real(0.0);
}

inline OpenHurricane::real OpenHurricane::transport::Di(const real p, const real T, const realArray &xi,
                                                const PtrList<transport> &tranPtr,
                                                realArray &Dij) const {
    return real(0.0);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::transport::Di(const real p, const real T,
                                                              const real Tv1p5, const realArray &xi,
                                                              const PtrList<transport> &tranPtr,
                                                              realArray &Dij) const {
    return real(0);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::transport::Di(const real p, const real T, const real Tv1p5, const realArray &xi,
                         const PtrList<transport> &tranPtr) const {
    return real(0);
}
