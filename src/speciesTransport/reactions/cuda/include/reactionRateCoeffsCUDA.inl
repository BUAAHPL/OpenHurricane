#include "reactionRateCoeffsCUDA.hpp"
/*!
 * \file reactionRateCoeffsCUDA.inl
 * \brief The In-Line functions of the <i>reactionRateCoeffsCUDA.hpp</i> file.
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

inline cu_host
OpenHurricane::cuChem::reactionCoeffs::reactionCoeffs()
    : size_(), index_(), stoichCoeff_(), order_() {}

inline cu_host
OpenHurricane::cuChem::reactionCoeffs::reactionCoeffs(
    const cu_ushort nrc, const cu_ushort maxSize,
    const cu_ushort *__restrict__ sizePtr,
    const cu_ushort *__restrict__ indexPtr,
    const cu_real *__restrict__ stoichCoeffPtr,
    const cu_real *__restrict__ orderPtr)
    : size_(nrc, sizePtr), index_(nrc, maxSize, indexPtr),
      stoichCoeff_(nrc, maxSize, stoichCoeffPtr),
      order_(nrc, maxSize, orderPtr) {}

inline cu_host
OpenHurricane::cuChem::reactionCoeffs::reactionCoeffs(
    const cu_ushort nrc, const cu_ushort maxSize,
    const cu_ushort *__restrict__ sizePtr,
    const cu_ushort *__restrict__ indexPtr,
    const cu_real *__restrict__ stoichCoeffPtr,
    const cu_real *__restrict__ orderPtr, const cudaStreams &streams)
    : size_(nrc, sizePtr, streams()), index_(nrc, maxSize, indexPtr),
      stoichCoeff_(nrc, maxSize, stoichCoeffPtr),
      order_(nrc, maxSize, orderPtr) {
    index_.copyFromHostAsync(indexPtr, streams());
    stoichCoeff_.copyFromHostAsync(stoichCoeffPtr, streams());
    order_.copyFromHostAsync(orderPtr, streams());
}

inline cu_host
OpenHurricane::cuChem::reactionCoeffs::reactionCoeffs(
    const cu_ushort nrc,
    const cuChemInterface::reactionCoeffsInterface &reacInt)
    : size_(nrc, reacInt.sizePtr_),
      index_(nrc, reacInt.maxSpcSizeInRec_, reacInt.reactionCoeffsIndexPtr_),
      stoichCoeff_(nrc, reacInt.maxSpcSizeInRec_, reacInt.stoichCoeffPtr_),
      order_(nrc, reacInt.maxSpcSizeInRec_, reacInt.orderPtr_) {}

inline cu_host
OpenHurricane::cuChem::reactionCoeffs::reactionCoeffs(
    const cu_ushort nrc,
    const cuChemInterface::reactionCoeffsInterface &reacInt,
    const cudaStreams &streams)
    : size_(nrc, reacInt.sizePtr_, streams()),
      index_(nrc, reacInt.maxSpcSizeInRec_),
      stoichCoeff_(nrc, reacInt.maxSpcSizeInRec_),
      order_(nrc, reacInt.maxSpcSizeInRec_) {
    index_.copyFromHostAsync(reacInt.reactionCoeffsIndexPtr_, streams());
    stoichCoeff_.copyFromHostAsync(reacInt.stoichCoeffPtr_, streams());
    order_.copyFromHostAsync(reacInt.orderPtr_, streams());
}

inline cu_dual
OpenHurricane::cuChem::reactionCoeffs::reactionCoeffs(
    const reactionCoeffs &rc)
    : size_(rc.size_), index_(rc.index_), stoichCoeff_(rc.stoichCoeff_),
      order_(rc.order_) {}

inline cu_dual
    OpenHurricane::cuChem::reactionCoeffs::~reactionCoeffs() noexcept {}

inline cu_host void
OpenHurricane::cuChem::reactionCoeffs::destroy() {
    destroyCuArray(size_);
    destroyCuArray(index_);
    destroyCuArray(stoichCoeff_);
    destroyCuArray(order_);
}

inline cu_device cu_ushort
OpenHurricane::cuChem::reactionCoeffs::size(
    const cu_ushort irc) const {
    return size_(irc);
}

inline cu_dual cu_real
OpenHurricane::cuChem::reactionRateCoeffs::ArrheniusReactionRate(
    const cu_real *__restrict__ a, const cu_real *__restrict__ b,
    const cu_real *__restrict__ Ta, const cu_integer ri,
    const cu_real T) const {
    return a[ri] * std::pow(T, b[ri]) * std::exp(-Ta[ri] / T);
}

inline OpenHurricane::cuChem::reactionRateCoeffs::reactionRateCoeffs()
    : a_(), b_(), Ta_(), destroyed_(false) {}

inline cu_host
OpenHurricane::cuChem::reactionRateCoeffs::reactionRateCoeffs(
    const cu_ushort nrc, const cu_real *__restrict__ aPtr,
    const cu_real *__restrict__ bPtr, const cu_real *__restrict__ TaPtr)
    : a_(nrc, aPtr), b_(nrc, bPtr), Ta_(nrc, TaPtr), destroyed_(false) {}

inline OpenHurricane::cuChem::reactionRateCoeffs::reactionRateCoeffs(
    const cu_ushort nrc, const cu_real *__restrict__ aPtr,
    const cu_real *__restrict__ bPtr, const cu_real *__restrict__ TaPtr,
    const cudaStreams &streams)
    : a_(nrc, aPtr, streams()), b_(nrc, bPtr, streams()),
      Ta_(nrc, TaPtr, streams()), destroyed_(false) {}

inline OpenHurricane::cuChem::reactionRateCoeffs::reactionRateCoeffs(
    const cu_ushort nrc,
    const cuChemInterface::reactionRateCoeffsInterface &reacInt)
    : a_(nrc, reacInt.aPtr_), b_(nrc, reacInt.bPtr_), Ta_(nrc, reacInt.TaPtr_),
      destroyed_(false) {}

inline OpenHurricane::cuChem::reactionRateCoeffs::reactionRateCoeffs(
    const cu_ushort nrc,
    const cuChemInterface::reactionRateCoeffsInterface &reacInt,
    const cudaStreams &streams)
    : a_(nrc, reacInt.aPtr_, streams()), b_(nrc, reacInt.bPtr_, streams()),
      Ta_(nrc, reacInt.TaPtr_, streams()), destroyed_(false) {}

inline cu_dual
OpenHurricane::cuChem::reactionRateCoeffs::reactionRateCoeffs(
    const reactionRateCoeffs &rrf)
    : a_(rrf.a_), b_(rrf.b_), Ta_(rrf.Ta_), destroyed_(false) {}

cu_dual inline OpenHurricane::cuChem::reactionRateCoeffs::
    ~reactionRateCoeffs() noexcept {}

inline cu_host void
OpenHurricane::cuChem::reactionRateCoeffs::destroy() {
    if (!destroyed_) {
        destroyCuArray(a_);
        destroyCuArray(b_);
        destroyCuArray(Ta_);
        destroyed_ = true;
    }
}

inline cu_device cu_real
OpenHurricane::cuChem::reactionRateCoeffs::k(const cu_ushort ri,
                                                    const cu_real T) const {
    return ArrheniusReactionRate(a_.dDataPtr(), b_.dDataPtr(), Ta_.dDataPtr(),
                                 ri, T);
}

#endif // CUDA_PARALLEL