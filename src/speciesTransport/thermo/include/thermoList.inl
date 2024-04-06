#include "thermoList.hpp"
/*!
 * \file thermoList.inl
 * \brief The In-Line functions of thermo table properties.
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

hur_nodiscard inline const OpenHurricane::PtrList<OpenHurricane::thermo> &
OpenHurricane::thermoList::thTable() const noexcept {
    return thList_;
}

hur_nodiscard inline const OpenHurricane::equationOfState &
OpenHurricane::thermoList::eos() const noexcept {
    return *eosPtr_;
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::p(const real rho, const real T,
                              const realArray &yi) const {
    return eosPtr_->pm(rho, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::p(const real rho, const real T,
                              const PtrList<cellRealArray> &yi,
                              const integer cellI) const {
    return eosPtr_->pm(rho, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::rho(const real p, const real T,
                                const realArray &yi) const {
    return eosPtr_->rhom(p, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::rho(const real p, const real T,
                                const PtrList<cellRealArray> &yi,
                                const integer cellI) const {
    return eosPtr_->rhom(p, T, yi, cellI);
}
hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cp0i(const real T, const integer i) const {
    return thList_[i].cp0(T);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::Cp0i(const real T, const integer i) const {
    return thList_[i].Cp0(T);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cp_rho(const real rho, const real T,
                                   const realArray &yi) const {
    return cp0(T, yi) + eosPtr_->cpm_rho(rho, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cp_rho(const real rho, const real T,
                                   const PtrList<cellRealArray> &yi,
                                   const integer cellI) const {
    return cp0(T, yi, cellI) + eosPtr_->cpm_rho(rho, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cp_p(const real p, const real T,
                                 const realArray &yi) const {
    return cp0(T, yi) + eosPtr_->cpm_p(p, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cp_p(const real p, const real T,
                                 const PtrList<cellRealArray> &yi,
                                 const integer cellI) const {
    return cp0(T, yi, cellI) + eosPtr_->cpm_p(p, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cpi_rho(const real rho, const real T,
                                    const integer i) const {
    return cp0i(T, i) + eosPtr_->cpi_rho(rho, T, i);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cpi_p(const real p, const real T,
                                  const integer i) const {
    return cp0i(T, i) + eosPtr_->cpi_p(p, T, i);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::Cp_rho(const real rho, const real T,
                                   const realArray &xi) const {
    return Cp0(T, xi) + eosPtr_->Cpm_rho(rho, T, xi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::Cp_p(const real p, const real T,
                                 const realArray &xi) const {
    return Cp0(T, xi) + eosPtr_->Cpm_p(p, T, xi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::Cpi_rho(const real rho, const real T,
                                    const integer i) const {
    return Cp0i(T, i) + eosPtr_->Cpi_rho(rho, T, i);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::Cpi_p(const real p, const real T,
                                  const integer i) const {
    return Cp0i(T, i) + eosPtr_->Cpi_p(p, T, i);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::Cv0i(const real T, const integer i) const {
    return thList_[i].Cv0(T);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cv_p(const real p, const real T,
                                 const realArray &yi) const {
    return cp_p(p, T, yi) - eosPtr_->cpMcvm_p(p, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cv_p(const real p, const real T,
                                 const PtrList<cellRealArray> &yi,
                                 const integer cellI) const {
    return cp_p(p, T, yi, cellI) - eosPtr_->cpMcvm_p(p, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cv_rho(const real rho, const real T,
                                   const realArray &yi) const {
    real p = eosPtr_->pm(rho, T, yi);
    return cv_p(p, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::cv_rho(const real rho, const real T,
                                   const PtrList<cellRealArray> &yi,
                                   const integer cellI) const {
    real p = eosPtr_->pm(rho, T, yi, cellI);
    return cv_p(p, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::gamma(const real p, const real T,
                                  const realArray &yi) const {
    const real cpt = cp_p(p, T, yi);
    return cpt / (cpt - eosPtr_->cpMcvm_p(p, T, yi));
    //return cp_p(p, T, yi) / cv_p(p, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::gamma(const real p, const real T,
                                  const PtrList<cellRealArray> &yi,
                                  const integer cellI) const {
    const real cpt = cp_p(p, T, yi, cellI);
    return cpt / (cpt - eosPtr_->cpMcvm_p(p, T, yi, cellI));
    //return cp_p(p, T, yi, cellI) / cv_p(p, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::gammai(const real p, const real T,
                                   const integer i) const {
    return thList_[i].gamma_p(p, T);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::gammaTh(const real p, const real rho, const real T,
                                    const realArray &yi) const {
    return gamma(p, T, yi) * eosPtr_->gammaThCorrectM(p, rho, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::gammaTh(const real p, const real rho, const real T,
                                    const PtrList<cellRealArray> &yi,
                                    const integer cellI) const {
    return gamma(p, T, yi, cellI) *
           eosPtr_->gammaThCorrectM(p, rho, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::gammaThi(const real p, const real rho, const real T,
                                     const integer i) const {
    return thList_[i].gammaThCorrect(rho, T);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::ha_p(const real p, const real T,
                                 const realArray &yi) const {
    return ha0(T, yi) + eosPtr_->hm_p(p, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::ha_p(const real p, const real T,
                                 const PtrList<cellRealArray> &yi,
                                 const integer cellI) const {
    return ha0(T, yi, cellI) + eosPtr_->hm_p(p, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::ha_p(const real p, const real T,
                                 const integer i) const {
    return thList_[i].ha0(T) + eosPtr_->hi_p(p, T, i);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::ha_rho(const real rho, const real T,
                                   const realArray &yi) const {
    return ha0(T, yi) + eosPtr_->hm_rho(rho, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::ha_rho(const real rho, const real T,
                                   const PtrList<cellRealArray> &yi,
                                   const integer cellI) const {
    return ha0(T, yi, cellI) + eosPtr_->hm_rho(rho, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::s_rho(const real rhom, const real T,
                                  const realArray &yi) const {
    real pm = eosPtr_->pm(rhom, T, yi);
    return s_p(pm, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::s_rho(const real rhom, const real T,
                                  const PtrList<cellRealArray> &yi,
                                  const integer cellI) const {
    real pm = eosPtr_->pm(rhom, T, yi, cellI);
    return s_p(pm, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::hs_p(const real p, const real T,
                                 const realArray &yi) const {
    return ha_p(p, T, yi) - hc(yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::hs_p(const real p, const real T,
                                 const PtrList<cellRealArray> &yi,
                                 const integer cellI) const {
    return ha_p(p, T, yi, cellI) - hc(yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::hs_rho(const real rho, const real T,
                                   const realArray &yi) const {
    return ha_rho(rho, T, yi) - hc(yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::hs_rho(const real rho, const real T,
                                   const PtrList<cellRealArray> &yi,
                                   const integer cellI) const {
    return ha_rho(rho, T, yi, cellI) - hc(yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::ea_p(const real p, const real T,
                                 const realArray &yi) const {
    return ha_p(p, T, yi) - p / this->rho(p, T, yi);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::ea_p(const real p, const real T,
                                 const PtrList<cellRealArray> &yi,
                                 const integer cellI) const {
    return ha_p(p, T, yi, cellI) - p / this->rho(p, T, yi, cellI);
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::ea_rho(const real rho, const real T,
                                   const realArray &yi) const {
    return ha_rho(rho, T, yi) - this->p(rho, T, yi) / rho;
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::thermoList::ea_rho(const real rho, const real T,
                                   const PtrList<cellRealArray> &yi,
                                   const integer cellI) const {
    return ha_rho(rho, T, yi, cellI) - this->p(rho, T, yi, cellI) / rho;
}

inline OpenHurricane::real
OpenHurricane::thermoList::THa_rho(real ha, real rho, real T0, integer &iFlag,
                                    const realArray &yi) const {
    return T(ha, rho, T0, iFlag, yi, &thermoList::ha_rho, &thermoList::cp_rho,
             &thermoList::limit);
}

inline OpenHurricane::real
OpenHurricane::thermoList::THa_rho(real ha, real rho, real T0, integer &iFlag,
                                    const PtrList<cellRealArray> &yi,
                                    const integer cellI) const {
    return T(ha, rho, T0, iFlag, yi, cellI, &thermoList::ha_rho,
             &thermoList::cp_rho, &thermoList::limit);
}

inline OpenHurricane::real
OpenHurricane::thermoList::THa_p(real ha, real p, real T0, integer &iFlag,
                                  const realArray &yi) const {
    return T(ha, p, T0, iFlag, yi, &thermoList::ha_p, &thermoList::cp_p,
             &thermoList::limit);
}
inline OpenHurricane::real
OpenHurricane::thermoList::THa_p(real ha, real p, real T0, integer &iFlag,
                                  const PtrList<cellRealArray> &yi,
                                  const integer cellI) const {
    return T(ha, p, T0, iFlag, yi, cellI, &thermoList::ha_p,
             &thermoList::cp_p, &thermoList::limit);
}
inline OpenHurricane::real
OpenHurricane::thermoList::TEa_rho(real ea, real rho, real T0, integer &iFlag,
                                    const realArray &yi) const {
    return T(ea, rho, T0, iFlag, yi, &thermoList::ea_rho, &thermoList::cv_rho,
             &thermoList::limit);
}

inline OpenHurricane::real
OpenHurricane::thermoList::TEa_rho(real ea, real rho, real T0, integer &iFlag,
                                    const PtrList<cellRealArray> &yi,
                                    const integer cellI) const {
    return T(ea, rho, T0, iFlag, yi, cellI, &thermoList::ea_rho,
             &thermoList::cv_rho, &thermoList::limit);
}

inline OpenHurricane::real
OpenHurricane::thermoList::TEa_p(real ea, real p, real T0, integer &iFlag,
                                  const realArray &yi) const {
    return T(ea, p, T0, iFlag, yi, &thermoList::ea_p, &thermoList::cv_p,
             &thermoList::limit);
}

inline OpenHurricane::real
OpenHurricane::thermoList::TEa_p(real ea, real p, real T0, integer &iFlag,
                                  const PtrList<cellRealArray> &yi,
                                  const integer cellI) const {
    return T(ea, p, T0, iFlag, yi, cellI, &thermoList::ea_p,
             &thermoList::cv_p, &thermoList::limit);
}
