#include "reactionRatePressureDependent.hpp"
/*!
 * \file reactionRatePressureDependent.inl
 * \brief The In-Line functions of the <i>reactionRatePressureDependent.hpp</i> file.
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

constexpr inline cu_dual cu_real
OpenHurricane::cuChem::pressureDependentCoefficients::fallOffFunctions::
    Lindemann(const cu_real T, const cu_real Pr) {
    return cu_real(1);
}

inline OpenHurricane::cuChem::pressureDependentCoefficients::
    pressureDependentCoefficients()
    : index_(), a_(), b_(), Ta_(), fallOffType_(), fallOffCoeff_() {}

inline OpenHurricane::cuChem::pressureDependentCoefficients::
    pressureDependentCoefficients(
        const cu_ushort nrc, const cu_ushort nrcPD,
        const cu_short *__restrict__ indexPtr,
        const cu_real *__restrict__ aPtr, const cu_real *__restrict__ bPtr,
        const cu_real *__restrict__ TaPtr,
        const cu_ushort *__restrict__ fallOffTypePtr,
        const cu_real *__restrict__ fallOffCoeffPtr)
    : index_(nrc, indexPtr), a_(nrcPD, aPtr), b_(nrcPD, bPtr),
      Ta_(nrcPD, TaPtr), fallOffType_(nrcPD, fallOffTypePtr),
      fallOffCoeff_(nrcPD, 5, fallOffCoeffPtr) {}

inline OpenHurricane::cuChem::pressureDependentCoefficients::
    pressureDependentCoefficients(
        const cu_ushort nrc, const cu_ushort nrcPD,
        const cu_short *__restrict__ indexPtr,
        const cu_real *__restrict__ aPtr, const cu_real *__restrict__ bPtr,
        const cu_real *__restrict__ TaPtr,
        const cu_ushort *__restrict__ fallOffTypePtr,
        const cu_real *__restrict__ fallOffCoeffPtr,
        const cudaStreams &streams)
    : index_(nrc, indexPtr, streams()), a_(nrcPD, aPtr, streams()),
      b_(nrcPD, bPtr, streams()), Ta_(nrcPD, TaPtr, streams()),
      fallOffType_(nrcPD, fallOffTypePtr, streams()), fallOffCoeff_(nrcPD, 5) {
    fallOffCoeff_.copyFromHostAsync(fallOffCoeffPtr, streams());
}

inline OpenHurricane::cuChem::pressureDependentCoefficients::
    pressureDependentCoefficients(
        const cu_ushort nrc,
        const cuChemInterface::pressureDepCoeffInterface &pressInt)
    : index_(nrc, pressInt.indexPtr_), a_(pressInt.nrcPD_, pressInt.aPtr_),
      b_(pressInt.nrcPD_, pressInt.bPtr_),
      Ta_(pressInt.nrcPD_, pressInt.TaPtr_),
      fallOffType_(pressInt.nrcPD_, pressInt.fallOffTypePtr_),
      fallOffCoeff_(pressInt.nrcPD_, 5, pressInt.fallOffCoeffPtr_) {}

inline OpenHurricane::cuChem::pressureDependentCoefficients::
    pressureDependentCoefficients(
        const cu_ushort nrc,
        const cuChemInterface::pressureDepCoeffInterface &pressInt,
        const cudaStreams &streams)
    : index_(nrc, pressInt.indexPtr_, streams()),
      a_(pressInt.nrcPD_, pressInt.aPtr_, streams()),
      b_(pressInt.nrcPD_, pressInt.bPtr_, streams()),
      Ta_(pressInt.nrcPD_, pressInt.TaPtr_, streams()),
      fallOffType_(pressInt.nrcPD_, pressInt.fallOffTypePtr_, streams()),
      fallOffCoeff_(pressInt.nrcPD_, 5) {
    fallOffCoeff_.copyFromHostAsync(pressInt.fallOffCoeffPtr_, streams());
}

inline cu_dual
OpenHurricane::cuChem::pressureDependentCoefficients::
    pressureDependentCoefficients(const pressureDependentCoefficients &pdc)
    : index_(pdc.index_), a_(pdc.a_), b_(pdc.b_), Ta_(pdc.Ta_),
      fallOffType_(pdc.fallOffType_), fallOffCoeff_(pdc.fallOffCoeff_) {}

cu_dual inline OpenHurricane::cuChem::
    pressureDependentCoefficients::~pressureDependentCoefficients() noexcept {}

inline cu_host void
OpenHurricane::cuChem::pressureDependentCoefficients::destroy() {
    destroyCuArray(index_);
    destroyCuArray(a_);
    destroyCuArray(b_);
    destroyCuArray(Ta_);
    destroyCuArray(fallOffType_);
    destroyCuArray(fallOffCoeff_);
}

inline cu_device cu_real
OpenHurricane::cuChem::pressureDependentCoefficients::k(
    const cu_ushort ri, const cu_real T) const {
    const auto rii = index_(ri);
    return a_(rii) * std::pow(T, b_(rii)) * std::exp(-Ta_(rii) / T);
}

inline cu_device cu_real
OpenHurricane::cuChem::pressureDependentCoefficients::Pr(
    const cu_real k0, const cu_real kinf, const cu_real M) const {
    return k0 / kinf * M;
}

inline cu_dual cu_real OpenHurricane::cuChem::
    pressureDependentCoefficients::UnimolecularFallOffRate(
        const cu_real kinf, const cu_real Pr, const cu_real F) const {
    return kinf * (Pr / (cu_real(1.0) + Pr)) * F;
}

inline cu_dual cu_real
OpenHurricane::cuChem::pressureDependentCoefficients::
    UnimolecularFallOffFactor(const cu_real Pr, const cu_real F) const {
    return (Pr / (cu_real(1.0) + Pr)) * F;
}

inline cu_dual cu_real OpenHurricane::cuChem::
    pressureDependentCoefficients::BimolecularFallOffRate(
        const cu_real k0, const cu_real Pr, const cu_real F) const {
    return k0 * (cu_real(1.0) / (cu_real(1.0) + Pr)) * F;
}

inline cu_dual cu_real
OpenHurricane::cuChem::pressureDependentCoefficients::
    BimolecularFallOffFactor(const cu_real Pr, const cu_real F) const {
    return (cu_real(1.0) / (cu_real(1.0) + Pr)) * F;
}

#endif // CUDA_PARALLEL
