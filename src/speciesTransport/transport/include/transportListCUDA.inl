#include "transportListCUDA.hpp"
/*!
 * \file transportListCUDA.inl
 * \brief The In-Line functions of the <i>transportListCUDA.hpp</i> file.
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

#ifdef CUDA_PARALLEL

inline cu_dual
OpenHurricane::CUDATransport::transportListCUDA::transportListCUDA(const transportListCUDA &tr)
    : transport_(tr.transport_) {}

inline cu_dual OpenHurricane::CUDATransport::transportListCUDA::~transportListCUDA() noexcept {
}

inline cu_host void OpenHurricane::CUDATransport::transportListCUDA::destroy() {
    transport_.destroy();
}

inline cu_dual const OpenHurricane::CUDATransport::kineticTheoryCUDA &
OpenHurricane::CUDATransport::transportListCUDA::transport() const noexcept {
    return transport_;
}

inline cu_device void OpenHurricane::CUDATransport::transportListCUDA::mu(
    cu_real *__restrict__ mui, const cu_real *__restrict__ xi,
    const cu_real *__restrict__ wci, const cu_ushort i) const {
    mui[i] *= wci[i];
}

inline cu_device void OpenHurricane::CUDATransport::transportListCUDA::muKappa(
    cu_real *__restrict__ mui, cu_real *__restrict__ kappai, const cu_real *__restrict__ xi,
    const cu_real *__restrict__ wci, const cu_ushort i) const {
    mui[i] *= wci[i];
    //kappai[i] *= (wci[i] / cu_real(1.065));
    kappai[i] *= wci[i];
}

inline cu_device cu_real OpenHurricane::CUDATransport::transportListCUDA::PhiMuij(
    const cu_real mui, const cu_real muj, const cu_real Wi, const cu_real Wj) const {
    const cu_real oneInv = cu_real(1) / sqrt(cu_real(8));
    return oneInv / sqrt(cu_real(1) + Wi / Wj) *
           cuSqr(cu_real(1) + sqrt(mui / muj * sqrt(Wj / Wi)));
}

inline cu_device cu_real OpenHurricane::CUDATransport::transportListCUDA::PhiKij(
    const cu_real mui, const cu_real muj, const cu_real Wi, const cu_real Wj) const {
    //return cu_real(1.065) * PhiMuij(mui, muj, Wi, Wj);
    return PhiMuij(mui, muj, Wi, Wj);
}

inline cu_device int OpenHurricane::CUDATransport::transportListCUDA::nsp() const {
    return transport_.nsp();
}

inline cu_device const OpenHurricane::cuChem::speciesTable &
OpenHurricane::CUDATransport::transportListCUDA::species() const noexcept {
    return transport_.species();
}

#endif // CUDA_PARALLEL