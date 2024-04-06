/*!
 * \file kineticTheory.cpp
 * \brief Main subroutines for transport properties by kinetic theory.
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

#include "kineticTheory.hpp"

namespace OpenHurricane {
    createClassNameStr(kineticTheory, "kinetic");
    registerObjFty(transport, kineticTheory, controller);
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::kineticTheory::mu(const real p, const real T) const {
    real Tstar = T * rekb_;
    const speciesList &sp = species();
    integer ind = index();

    return muPre_ * sqrt(sp[ind].W() * T) / (sigma2_ * omega22(Tstar));
}

hur_nodiscard OpenHurricane::real OpenHurricane::kineticTheory::kappa(const real p, const real T,
                                                              const real cpi) const {
    const speciesList &sp = species();
    integer ind = index();
    return mu(p, T) * sp[ind].Ri() * (cpi / sp[ind].Ri() + 1.25);
}

hur_nodiscard OpenHurricane::real
OpenHurricane::kineticTheory::kappa(const real p, const real T, const real mui, const real cpi) const {
    const speciesList &sp = species();
    integer ind = index();
    return mui * sp[ind].Ri() * (cpi / sp[ind].Ri() + 1.25);
}

hur_nodiscard OpenHurricane::real OpenHurricane::kineticTheory::Di(const real p, const real T,
                                                           const real Tv1p5, const realArray &xi,
                                                           const PtrList<transport> &tranPtr,
                                                           realArray &Dij) const {
    const speciesList &sp = species();
    const integer ind = index();
    const integer nsp = xi.size();

    const auto pp = p;
    Dij = Zero;
    if (DijConstPartPtr_ == nullptr) {
        const real invW = real(1.0) / sp[ind].W();
        DijConstPartPtr_ = new realArray(nsp);
        for (integer j = 0; j < nsp; ++j) {
            const kineticTheory *ktj = static_cast<const kineticTheory *>(tranPtr(j));
            const real sigmaij = 0.5 * (sigma_ + ktj->sigma_);

            (*DijConstPartPtr_)[j] = 1.858e-7 * sqrt(invW + real(1.0) / sp[j].W()) / sqr(sigmaij) *
                                     constant::physicalConstant::Patm;
        }
    }

    if (rekbij05Ptr_ == nullptr) {
        rekbij05Ptr_ = new realArray(nsp);
        for (integer j = 0; j < nsp; ++j) {
            const kineticTheory *ktj = static_cast<const kineticTheory *>(tranPtr(j));
            (*rekbij05Ptr_)[j] = inv(sqrt(ekb_ * ktj->ekb_));
        }
    }

    for (integer j = 0; j < nsp; ++j) {
        if (ind != j) {
            const kineticTheory *ktj = static_cast<const kineticTheory *>(tranPtr(j));
            //const real Tstar = T / sqrt(ekb_ * ktj->ekb_);
            const real Tstar = T * (*rekbij05Ptr_)[j];
            const real Tss = 1.0 / (sqr(Tstar + real(0.5)));
            const real OmegaDij = exp(-0.145 * log(Tstar)) + Tss;

            Dij[j] = Tv1p5 * (*DijConstPartPtr_)[j] / (pp * OmegaDij);
        }
    }
    real sum = Zero;

    for (integer j = 0; j < nsp; ++j) {
        if (ind != j) {
            sum += xi[j] / Dij[j];
        }
    }

    real Dim = Zero;
    if (xi[ind] < 1.0 - veryTiny) {
        Dim = (1.0 - xi[ind]) / sum;
    }
    if (isinf(Dim)) {
        return 0;
    }
    return Dim;
}

hur_nodiscard OpenHurricane::real
OpenHurricane::kineticTheory::Di(const real p, const real T, const real Tv1p5, const realArray &xi,
                             const PtrList<transport> &tranPtr) const {
    const speciesList &sp = species();
    const integer ind = index();
    const integer nsp = xi.size();

    const auto pp = p;
    if (DijConstPartPtr_ == nullptr) {
        const real invW = real(1.0) / sp[ind].W();
        DijConstPartPtr_ = new realArray(nsp);
        for (integer j = 0; j < nsp; ++j) {
            const kineticTheory *ktj = static_cast<const kineticTheory *>(tranPtr(j));
            const real sigmaij = 0.5 * (sigma_ + ktj->sigma_);

            (*DijConstPartPtr_)[j] = 1.858e-7 * sqrt(invW + real(1.0) / sp[j].W()) / sqr(sigmaij) *
                                     constant::physicalConstant::Patm;
        }
    }

    if (rekbij05Ptr_ == nullptr) {
        rekbij05Ptr_ = new realArray(nsp);
        for (integer j = 0; j < nsp; ++j) {
            const kineticTheory *ktj = static_cast<const kineticTheory *>(tranPtr(j));
            (*rekbij05Ptr_)[j] = inv(sqrt(ekb_ * ktj->ekb_));
        }
    }

    real sum = Zero;
    for (integer j = 0; j < nsp; ++j) {
        if (ind != j) {
            const kineticTheory *ktj = static_cast<const kineticTheory *>(tranPtr(j));
            //const real Tstar = T / sqrt(ekb_ * ktj->ekb_);
            const real Tstar = T * (*rekbij05Ptr_)[j];
            const real Tss = 1.0 / (sqr(Tstar + real(0.5)));
            const real OmegaDij = exp(-0.145 * log(Tstar)) + Tss;

            real Dij = Tv1p5 * (*DijConstPartPtr_)[j] / (pp * OmegaDij);
            sum += xi[j] / Dij;
        }
    }

    real Dim = Zero;
    if (xi[ind] < 1.0 - veryTiny) {
        Dim = (1.0 - xi[ind]) / sum;
    }
    if (isinf(Dim)) {
        return 0;
    }
    return Dim;
}

OpenHurricane::kineticTheory &OpenHurricane::kineticTheory::operator=(const kineticTheory &kT) {
    if (this != std::addressof(kT)) {
        transport::operator=(kT);

        geometryIndex_ = kT.geometryIndex_;
        ekb_ = kT.ekb_;
        rekb_ = kT.rekb_;
        sigma_ = kT.sigma_;
        sigma2_ = kT.sigma2_;
        mmu_ = kT.mmu_;
        alpha_ = kT.alpha_;
        Zrot_ = kT.Zrot_;
    }
    return *this;
}
