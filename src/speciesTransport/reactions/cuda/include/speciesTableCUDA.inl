#include "speciesTableCUDA.hpp"
/*!
 * \file speciesTableCUDA.inl
 * \brief The In-Line functions of the <i>speciesTableCUDA.hpp</i> file.
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

inline cu_host OpenHurricane::cuChem::speciesTable::speciesTable()
    : W_(), thermo_() {}

inline cu_host OpenHurricane::cuChem::speciesTable::speciesTable(
    const cu_ushort nsp, const cu_real *__restrict__ WPtr,
    const cu_real *__restrict__ TCommonPtr,
    const cu_real *__restrict__ TLowPtr,
    const cu_real *__restrict__ THighPtr,
    const cu_real *__restrict__ highCpCoeffsPtr,
    const cu_real *__restrict__ lowCpCoeffsPtr)
    : W_(nsp, WPtr), thermo_(nsp, highCpCoeffsPtr, lowCpCoeffsPtr, TLowPtr,
                             THighPtr, TCommonPtr) {}

inline OpenHurricane::cuChem::speciesTable::speciesTable(
    const cu_ushort nsp, const cu_real *__restrict__ WPtr,
    const cu_real *__restrict__ TCommonPtr,
    const cu_real *__restrict__ TLowPtr,
    const cu_real *__restrict__ THighPtr,
    const cu_real *__restrict__ highCpCoeffsPtr,
    const cu_real *__restrict__ lowCpCoeffsPtr, const cudaStreams &streams)
    : W_(nsp, WPtr, streams()),
      thermo_(nsp, highCpCoeffsPtr, lowCpCoeffsPtr, TLowPtr, THighPtr,
              TCommonPtr, streams) {}

inline OpenHurricane::cuChem::speciesTable::speciesTable(
    const cu_ushort nsp,
    const cuChemInterface::speciesTableInterface &sptInt)
    : W_(nsp, sptInt.WPtr_),
      thermo_(nsp, sptInt.highCpCoeffsPtr_, sptInt.lowCpCoeffsPtr_,
              sptInt.TLowPtr_, sptInt.THighPtr_, sptInt.TCommonPtr_) {}

inline OpenHurricane::cuChem::speciesTable::speciesTable(
    const cu_ushort nsp,
    const cuChemInterface::speciesTableInterface &sptInt,
    const cudaStreams &streams)
    : W_(nsp, sptInt.WPtr_, streams()),
      thermo_(nsp, sptInt.highCpCoeffsPtr_, sptInt.lowCpCoeffsPtr_,
              sptInt.TLowPtr_, sptInt.THighPtr_, sptInt.TCommonPtr_, streams) {}

inline cu_dual OpenHurricane::cuChem::speciesTable::speciesTable(
    const speciesTable &spt)
    : W_(spt.W_), thermo_(spt.thermo_) {}

inline cu_dual
    OpenHurricane::cuChem::speciesTable::~speciesTable() noexcept {}

cu_host void OpenHurricane::cuChem::speciesTable::destroy() {
    destroyCuArray(W_);
    thermo_.destroy();
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::nGdRT(const cu_ushort isp,
                                                  const cu_real Td) const {
    return thermo_.SH(Td, isp);
}

inline cu_dual void OpenHurricane::cuChem::speciesTable::yiToci(
    const cu_real rho, const cu_real *__restrict__ yi,
    const cu_real *__restrict__ wi, const cu_ushort isp,
    cu_real *__restrict__ ci) const {
    ci[isp] = rho * yi[isp] / wi[isp];
}

inline cu_device void OpenHurricane::cuChem::speciesTable::yiToci(
    const cu_real rho, const cu_real *__restrict__ yi,
    const cu_ushort isp, cu_real *__restrict__ ci) const {
    ci[isp] = rho * yi[isp] / W_(isp);
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::yiToci(
    const cu_real rho, const cu_real yi, const cu_ushort isp) const {
    return rho * yi / Wi(isp);
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::Wi(const cu_ushort isp) const {
    return W_(isp);
}

inline cu_device const OpenHurricane::CUDAThermo::JANAFCUDA &
OpenHurricane::cuChem::speciesTable::thermo() const {
    return thermo_;
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::cp0(const cu_real T,
                                                const cu_ushort isp) const {
    return thermo_.Cp0(T, isp) / Wi(isp);
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::Cp0(const cu_real T,
                                                const cu_ushort isp) const {
    return thermo_.Cp0(T, isp);
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::cv0(const cu_real T,
                                                const cu_ushort isp) const {
    return cp0(T, isp) - Ri(isp);
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::Cv0(const cu_real T,
                                                const cu_ushort isp) const {
    return Cp0(T, isp) - Rgen;
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::ha0(const cu_real T,
                                                const cu_ushort isp) const {
    return thermo_.Ha0(T, isp) / Wi(isp);
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::Ha0(const cu_real T,
                                                const cu_ushort isp) const {
    return thermo_.Ha0(T, isp);
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::ea0(const cu_real T,
                                                const cu_ushort isp) const {
    return ha0(T, isp) - Ri(isp) * T;
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::Ea0(const cu_real T,
                                                const cu_ushort isp) const {
    return Ha0(T, isp) - Rgen * T;
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::s0(const cu_real T,
                                               const cu_ushort isp) const {
    return thermo_.S0(T, isp) / Wi(isp);
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::S0(const cu_real T,
                                               const cu_ushort isp) const {
    return thermo_.S0(T, isp);
}

inline cu_device cu_real
OpenHurricane::cuChem::speciesTable::Ri(const cu_ushort isp) const {
    return Rgen / Wi(isp);
}

inline cu_device void OpenHurricane::cuChem::speciesTable::yidWi(
    const cu_ushort isp, const cu_real *__restrict__ yi,
    cu_real *__restrict__ ydwi) const {
    ydwi[isp] = yi[isp] / Wi(isp);
}

inline cu_device void OpenHurricane::cuChem::speciesTable::yi2xi(
    const cu_ushort isp, const cu_real *__restrict__ yi, const cu_real wm,
    cu_real *__restrict__ xi) const {
    xi[isp] = yi[isp] * wm / Wi(isp);
    xi[isp] = cu_min(1, xi[isp]);
}

#endif // CUDA_PARALLEL