#include "transportList.hpp"
/*!
 * \file transportList.inl
 * \brief The In-Line functions of transport properties table.
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

hur_nodiscard inline const OpenHurricane::PtrList<OpenHurricane::transport> &
OpenHurricane::transportList::tranTable() const noexcept {
    return tranList_;
}

inline void OpenHurricane::transportList::muf(const real p, const real T, realArray &mui) const {
    for (integer i = 0; i < species_.size(); ++i) {
        mui[i] = mu(p, T, i);
    }
}

inline void OpenHurricane::transportList::kappaf(const real p, const real T, const realArray &mui,
                                              const realArray &cpi, realArray &kappaf) const {
    for (integer i = 0; i < species_.size(); ++i) {
        kappaf[i] = kappa(p, T, mui[i], cpi[i], i);
    }
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::transportList::D(const real p, const real T,
                                                                  const realArray &yi,
                                                                  const integer i) const {
    return tranList_[i].Di(p, T, species_.Yi2Xi(yi), tranList_);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::transportList::D_Xi(const real p, const real T,
                                                                     const realArray &xi,
                                                                     const integer i) const {
    return tranList_[i].Di(p, T, xi, tranList_);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::transportList::D(const real p, const real T,
                                                                  const realArray &yi,
                                                                  const integer i,
                                                                  realArray &Dij) const {
    return tranList_[i].Di(p, T, species_.Yi2Xi(yi), tranList_, Dij);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::transportList::D_Xi(const real p, const real T,
                                                                     const realArray &xi,
                                                                     const integer i,
                                                                     realArray &Dij) const {
    return tranList_[i].Di(p, T, xi, tranList_, Dij);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::transportList::D_Xi(const real p, const real T, const real Tv1p5, const realArray &xi,
                                const integer i, realArray &Dij) const {
    return tranList_[i].Di(p, T, Tv1p5, xi, tranList_, Dij);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::transportList::D_Xi(const real p, const real T,
                                                                     const real Tv1p5,
                                                                     const realArray &xi,
                                                                     const integer i) const {
    return tranList_[i].Di(p, T, Tv1p5, xi, tranList_);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::transportList::D(const real p, const real T,
                                                                  const PtrList<cellRealArray> &yi,
                                                                  const integer cellI,
                                                                  const integer i) const {
    return tranList_[i].Di(p, T, species_.Yi2Xi(yi, cellI), tranList_);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::transportList::D(const real p, const real T, const PtrList<cellRealArray> &yi,
                             const integer cellI, const integer i, realArray &Dij) const {
    return tranList_[i].Di(p, T, species_.Yi2Xi(yi, cellI), tranList_, Dij);
}