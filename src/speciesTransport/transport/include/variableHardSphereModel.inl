#include "variableHardSphereModel.hpp"
/*!
 * \file variableHardSphereModel.inl
 * \brief The In-Line functions of transport properties by the variable hard sphere model.
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

inline OpenHurricane::variableHardSphereModel::variableHardSphereModel(const speciesList &sp,
                                                                   const integer index,
                                                                   const real omega,
                                                                   const real Tref, const real dref,
                                                                   const real Prl)
    : transport(sp, index, Prl), omega_(omega), Tref_(Tref), dref_(dref), m_(), muref_() {
    m_ = sp.W(index) / 1000.0 / constant::physicalConstant::NA;
    calcMuRef();
}

inline OpenHurricane::variableHardSphereModel::variableHardSphereModel(const speciesList &sp,
                                                                   const integer index,
                                                                   const controller &cont)
    : transport(sp, index, cont), omega_(cont.findOrDefault<real>("omega", 0.0)),
      Tref_(cont.findOrDefault<real>("Tref", 300.0)),
      dref_(cont.findOrDefault<real>("dref", 1e-12)), m_(), muref_() {
    m_ = sp.W(index) / 1000.0 / constant::physicalConstant::NA;
    calcMuRef();
}

inline OpenHurricane::variableHardSphereModel::variableHardSphereModel(
    const variableHardSphereModel &tra)
    : transport(tra), omega_(tra.omega_), Tref_(tra.Tref_), dref_(tra.dref_), m_(tra.m_),
      muref_(tra.muref_) {}

inline OpenHurricane::variableHardSphereModel::variableHardSphereModel(
    const variableHardSphereModel &tra, const speciesList &sp)
    : transport(tra, sp), omega_(tra.omega_), Tref_(tra.Tref_), dref_(tra.dref_), m_(tra.m_),
      muref_(tra.muref_) {}

inline OpenHurricane::variableHardSphereModel::~variableHardSphereModel() noexcept {}

inline OpenHurricane::real OpenHurricane::variableHardSphereModel::mu(const real p, const real T) const {
    return muref_ * pow(T / Tref_, omega_);
}

inline OpenHurricane::real OpenHurricane::variableHardSphereModel::kappa(const real p, const real T,
                                                                 const real cpi) const {
    return mu(p, T) * cpi / Pr();
}

inline OpenHurricane::real OpenHurricane::variableHardSphereModel::kappa(const real p, const real T,
                                                                 const real mui,
                                                                 const real cpi) const {
    return mui * cpi / Pr();
}
