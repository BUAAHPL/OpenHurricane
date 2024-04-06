#include "kineticTheory.hpp"
/*!
 * \file kineticTheory.inl
 * \brief The In-Line functions of transport properties by kinetic theory.
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

hur_nodiscard inline OpenHurricane::real OpenHurricane::kineticTheory::omega22(const real Tstar) const {
    return omegaPre_ * exp(omegaExponent_ * log(Tstar)) + real(1.0) / sqr(Tstar + real(0.5));
}

inline OpenHurricane::kineticTheory::kineticTheory(const speciesList &sp, const integer index,
                                               const real Prl, const integer gi, const real ekb,
                                               const real sig, const real mmu, const real alp,
                                               const real Zr)
    : transport(sp, index, Prl), geometryIndex_(gi), ekb_(ekb), sigma_(sig), mmu_(mmu), alpha_(alp),
      Zrot_(Zr), DijConstPartPtr_(nullptr), rekbij05Ptr_(nullptr) {
    rekb_ = 1.0 / ekb_;
    sigma2_ = sqr(sigma_);
}

inline OpenHurricane::kineticTheory::kineticTheory(const speciesList &sp, const integer index,
                                               const controller &cont)
    : transport(sp, index, cont), geometryIndex_(cont.findType<integer>("geonetryIndex", 0)),
      ekb_(cont.findType<real>("ekb", 0.0)), sigma_(cont.findType<real>("sigma", 0.0)),
      mmu_(cont.findType<real>("mmu", 0.0)), alpha_(cont.findType<real>("alpha", 0.0)),
      Zrot_(cont.findType<real>("Zrot", 0.0)), DijConstPartPtr_(nullptr), rekbij05Ptr_(nullptr) {
    rekb_ = 1.0 / ekb_;
    sigma2_ = sqr(sigma_);
}

inline OpenHurricane::kineticTheory::kineticTheory(const kineticTheory &tra)
    : transport(tra), geometryIndex_(tra.geometryIndex_), ekb_(tra.ekb_), rekb_(tra.rekb_),
      sigma_(tra.sigma_), sigma2_(tra.sigma2_), mmu_(tra.mmu_), alpha_(tra.alpha_),
      Zrot_(tra.Zrot_), DijConstPartPtr_(nullptr), rekbij05Ptr_(nullptr) {
    if (tra.DijConstPartPtr_ != nullptr) {
        DijConstPartPtr_ = new realArray(species().size());
        *DijConstPartPtr_ = *tra.DijConstPartPtr_;
    }
    if (tra.rekbij05Ptr_ != nullptr) {
        rekbij05Ptr_ = new realArray(species().size());
        *rekbij05Ptr_ = *tra.rekbij05Ptr_;
    }
}

inline OpenHurricane::kineticTheory::kineticTheory(const kineticTheory &tra, const speciesList &sp)
    : transport(tra, sp), geometryIndex_(tra.geometryIndex_), ekb_(tra.ekb_), rekb_(tra.rekb_),
      sigma_(tra.sigma_), sigma2_(tra.sigma2_), mmu_(tra.mmu_), alpha_(tra.alpha_),
      Zrot_(tra.Zrot_), DijConstPartPtr_(nullptr), rekbij05Ptr_(nullptr) {
    if (tra.DijConstPartPtr_ != nullptr) {
        DijConstPartPtr_ = new realArray(species().size());
        *DijConstPartPtr_ = *tra.DijConstPartPtr_;
    }
    if (tra.rekbij05Ptr_ != nullptr) {
        rekbij05Ptr_ = new realArray(species().size());
        *rekbij05Ptr_ = *tra.rekbij05Ptr_;
    }
}

inline OpenHurricane::kineticTheory::~kineticTheory() noexcept {
    HurDelete(DijConstPartPtr_);
    HurDelete(rekbij05Ptr_);
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::kineticTheory::geometryIndex() const noexcept {
    return geometryIndex_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::kineticTheory::ekb() const noexcept {
    return ekb_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::kineticTheory::sigma() const noexcept {
    return sigma_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::kineticTheory::mmu() const noexcept {
    return mmu_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::kineticTheory::alpha() const noexcept {
    return alpha_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::kineticTheory::Zrot() const noexcept {
    return Zrot_;
}

hur_nodiscard inline OpenHurricane::real
OpenHurricane::kineticTheory::Di(const real p, const real T, const realArray &xi,
                             const PtrList<transport> &tranPtr) const {
    realArray Dij(species().size());
    return Di(p, T, xi, tranPtr, Dij);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::kineticTheory::Di(const real p, const real T,
                                                                  const realArray &xi,
                                                                  const PtrList<transport> &tranPtr,
                                                                  realArray &Dij) const {
    return Di(p, T, T * sqrt(T), xi, tranPtr, Dij);
}
