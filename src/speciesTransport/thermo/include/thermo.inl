#include "thermo.hpp"
/*!
 * \file thermo.inl
 * \brief The In-Line functions of thermo properties.
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

inline OpenHurricane::thermo::~thermo() noexcept {
    eosPtr_ = nullptr;
}

inline const OpenHurricane::equationOfState &OpenHurricane::thermo::eos() const {
    return *eosPtr_;
}

inline OpenHurricane::real OpenHurricane::thermo::inteCp0dT(const real T1, const real T2) const {
    return s0(T2) - s0(T1);
}

inline OpenHurricane::real OpenHurricane::thermo::Cp0(const real T) const {
    return cp0(T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::cp_p(const real pi, const real T) const {
    return cp0(T) + eosPtr_->cpi_p(pi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::Dcp_pDT(const real pi, const real T) const {
    return real(0);
}

inline OpenHurricane::real OpenHurricane::thermo::Cp_p(const real pi, const real T) const {
    return Cp0(T) + eosPtr_->Cpi_p(pi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::DCp_pDT(const real pi, const real T) const {
    return this->Dcp_pDT(pi, T) * this->Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::cp_rho(const real rhoi, const real T) const {
    return cp0(T) + eosPtr_->cpi_rho(rhoi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::Cp_rho(const real rhoi, const real T) const {
    return Cp0(T) + eosPtr_->Cpi_rho(rhoi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::cv0(const real T) const {
    return cp0(T) - eosPtr_->species()[speciesIndex_].Ri();
}

inline OpenHurricane::real OpenHurricane::thermo::Cv0(const real T) const {
    return cv0(T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::cv_p(const real pi, const real T) const {
    return cp_p(pi, T) - eosPtr_->cpMcvi_p(pi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::Dcv_pDT(const real pi, const real T) const {
    return Dcp_pDT(pi, T);
}

inline OpenHurricane::real OpenHurricane::thermo::Cv_p(const real pi, const real T) const {
    return cv_p(pi, T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::DCv_pDT(const real pi, const real T) const {
    return Dcv_pDT(pi, T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::cv_rho(const real rhoi, const real T) const {
    return cp_rho(rhoi, T) - eosPtr_->cpMcvi_rho(rhoi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::gamma_p(const real pi, const real T) const {
    const real cp = cp_p(pi, T);
    const real cv = cp - eosPtr_->cpMcvi_p(pi, T, speciesIndex_);
    return cp / cv;
}

inline OpenHurricane::real OpenHurricane::thermo::gamma_rho(const real rhoi, const real T) const {
    const real cp = cp_rho(rhoi, T);
    const real cv = cp - eosPtr_->cpMcvi_rho(rhoi, T, speciesIndex_);
    return cp / cv;
}

inline OpenHurricane::real OpenHurricane::thermo::gammaThCorrect(const real rhoi, const real T) const {
    return gamma_rho(rhoi, T) * eosPtr_->gammaThCorrect(rhoi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::u0(const real T) const {
    return ha0(T) - Ri() * T;
}

inline OpenHurricane::real OpenHurricane::thermo::U0(const real T) const {
    return u0(T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::Ha0(const real T) const {
    return ha0(T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::DHa0DT(const real T) const {
    return Dha0DT(T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::ha_p(const real pi, const real T) const {
    return ha0(T) + eosPtr_->hi_p(pi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::ha_rho(const real rhoi, const real T) const {
    return ha0(T) + eosPtr_->hi_rho(rhoi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::Ha_p(const real pi, const real T) const {
    return ha_p(pi, T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::Ha_rho(const real rhoi, const real T) const {
    return ha_rho(rhoi, T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::hs_p(const real pi, const real T) const {
    return ha_p(pi, T) - hc();
}

inline OpenHurricane::real OpenHurricane::thermo::Hs_p(const real pi, const real T) const {
    return hs_p(pi, T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::hs_rho(const real rhoi, const real T) const {
    return ha_rho(rhoi, T) - hc();
}

inline OpenHurricane::real OpenHurricane::thermo::Hs_rho(const real rhoi, const real T) const {
    return hs_rho(rhoi, T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::Hc() const {
    return hc() * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::ea(const real rho, const real p, const real T) const {
    //return ha_rho(rho, T) - p / rho;
    return ha_rho(rho, T) - Ri() * T;
}

inline OpenHurricane::real OpenHurricane::thermo::S0(const real T) {
    return s0(T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::DS0DT(const real T) {
    return Ds0DT(T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::s_p(const real pi, const real T) const {
    return s0(T) + eosPtr_->si_p(pi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::s_rho(const real rhoi, const real T) const {
    return s0(T) + eosPtr_->si_rho(rhoi, T, speciesIndex_);
}

inline OpenHurricane::real OpenHurricane::thermo::g0(const real T) const {
    return ha0(T) - T * s0(T);
}

inline OpenHurricane::real OpenHurricane::thermo::G0(const real T) const {
    return g0(T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::Dg0DT(const real T) const {
    return Dha0DT(T) - (s0(T) + T * Ds0DT(T));
}

inline OpenHurricane::real OpenHurricane::thermo::DG0DT(const real T) const {
    return Dg0DT(T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::a0(const real T) const {
    return u0(T) - T * s0(T);
}

inline OpenHurricane::real OpenHurricane::thermo::A0(const real T) const {
    return a0(T) * Wi();
}

inline OpenHurricane::real OpenHurricane::thermo::Ri() const {
    return eosPtr_->species()[speciesIndex_].Ri();
}

inline OpenHurricane::real OpenHurricane::thermo::Wi() const {
    return eosPtr_->species()[speciesIndex_].W();
}
