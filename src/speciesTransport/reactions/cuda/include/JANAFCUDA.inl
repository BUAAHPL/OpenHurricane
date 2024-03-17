#include "JANAFCUDA.hpp"
/*!
 * \file JANAFCUDA.inl
 * \brief The In-Line functions of the <i>JANAFCUDA.hpp</i> file.
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

#ifdef CUDA_PARALLEL

inline cu_host OpenHurricane::CUDAThermo::JANAFCUDA::JANAFCUDA()
    : highCpCoeffs_(), lowCpCoeffs_(), TLow_(), THigh_(), TCommon_() {}

inline cu_dual OpenHurricane::CUDAThermo::JANAFCUDA::JANAFCUDA(const JANAFCUDA &jc)
    : highCpCoeffs_(jc.highCpCoeffs_), lowCpCoeffs_(jc.lowCpCoeffs_), TLow_(jc.TLow_),
      THigh_(jc.THigh_), TCommon_(jc.TCommon_) {}

inline cu_dual OpenHurricane::CUDAThermo::JANAFCUDA::~JANAFCUDA() noexcept {}

inline cu_host void OpenHurricane::CUDAThermo::JANAFCUDA::destroy() {
    destroyCuArray(highCpCoeffs_);
    destroyCuArray(lowCpCoeffs_);
    destroyCuArray(TLow_);
    destroyCuArray(THigh_);
    destroyCuArray(TCommon_);
}

inline cu_device cu_real
OpenHurricane::CUDAThermo::JANAFCUDA::S0(const cu_real T, const cu_ushort isp) const {
    const auto &a = coeff(T, isp);
    return Rgen *
           ((((a(4, isp) / 4.0 * T + a(3, isp) / 3.0) * T + a(2, isp) / 2.0) * T + a(1, isp)) * T +
            a(0, isp) * log(T) + a(6, isp));
}

inline cu_device cu_real
OpenHurricane::CUDAThermo::JANAFCUDA::SH(const cu_real T, const cu_ushort isp) const {
    const auto &a = coeff(T, isp);

    return a(0, isp) * (log(T) - cu_real(1)) +
           (((a(4, isp) / cu_real(20) * T + a(3, isp) / cu_real(12)) * T +
             a(2, isp) / cu_real(6)) *
                T +
            a(1, isp) / cu_real(2)) *
               T +
           a(6, isp) - a(5, isp) / cu_max(T, cu_tiny);
}

inline cu_device const OpenHurricane::cu2DArray<cu_real> &
OpenHurricane::CUDAThermo::JANAFCUDA::coeff(const cu_real T, const cu_ushort isp) const {
    if (T < TCommon_(isp)) {
        return lowCpCoeffs_;
    } else {
        return highCpCoeffs_;
    }
}

inline cu_device cu_real
OpenHurricane::CUDAThermo::JANAFCUDA::limit(const cu_real T, const cu_ushort isp) const {
    if (T < TLow_(isp)) {
        return TLow_(isp);
    } else if (T > THigh_(isp)) {
        return THigh_(isp);
    }
    return T;
}

#endif // CUDA_PARALLEL