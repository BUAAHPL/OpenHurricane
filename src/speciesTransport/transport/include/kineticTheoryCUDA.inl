#include "kineticTheoryCUDA.hpp"
/*!
 * \file reactionTableCUDA.inl
 * \brief The In-Line functions of the <i>reactionTableCUDA.hpp</i> file.
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
OpenHurricane::CUDATransport::kineticTheoryCUDA::kineticTheoryCUDA(const kineticTheoryCUDA &jc)
    : ekb_(jc.ekb_), sigma_(jc.sigma_), species_(jc.species_) {}

inline cu_dual OpenHurricane::CUDATransport::kineticTheoryCUDA::~kineticTheoryCUDA() noexcept {}

inline cu_host void OpenHurricane::CUDATransport::kineticTheoryCUDA::destroy() {
    destroyCuArray(ekb_);
    destroyCuArray(sigma_);
}

inline cu_device cu_real OpenHurricane::CUDATransport::kineticTheoryCUDA::kappa(
    const cu_real T, const cu_real mui, const cu_real cpi, const cu_ushort isp) const {
    const auto Ri = species_.Ri(isp);
    return mui * Ri * (cpi / Ri + cu_real(1.25));
}

inline cu_device cu_real OpenHurricane::CUDATransport::kineticTheoryCUDA::omega22(
    const cu_real Tstar, const cu_ushort isp) const {
    return cu_real(1.147) * exp(cu_real(-0.145) * log(Tstar)) +
           cu_real(1.0) / cuSqr(Tstar + cu_real(0.5));
}

inline cu_device cu_real OpenHurricane::CUDATransport::kineticTheoryCUDA::omegaDij(
    const cu_real T, const cu_real ekbi, const cu_real ekbj) const {
    const cu_real Tstar = T / sqrt(ekbi * ekbj);
    return exp(cu_real(-0.145) * log(Tstar)) + cu_real(1.0) / (cuSqr(Tstar + cu_real(0.5)));
}

inline cu_device int OpenHurricane::CUDATransport::kineticTheoryCUDA::nsp() const {
    return ekb_.size();
}

inline cu_device const OpenHurricane::cuChem::speciesTable &
OpenHurricane::CUDATransport::kineticTheoryCUDA::species() const noexcept {
    return species_;
}

#endif // CUDA_PARALLEL