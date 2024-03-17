#include "reactionRateThirdBody.hpp"
/*!
 * \file reactionRateThirdBody.inl
 * \brief The In-Line functions of the <i>reactionRateThirdBody.hpp</i> file.
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

inline OpenHurricane::cuChem::thirdBodyCoefficients::
    thirdBodyCoefficients()
    : index_(), coefThird_() {}

inline cu_host
OpenHurricane::cuChem::thirdBodyCoefficients::thirdBodyCoefficients(
    const cu_ushort nsp, const cu_ushort nrc, const cu_ushort nrcThird,
    const cu_short *__restrict__ indexPtr,
    const cu_real *__restrict__ coefThirdPtr)
    : index_(nrc, indexPtr), coefThird_(nrcThird, nsp, coefThirdPtr) {}

inline OpenHurricane::cuChem::thirdBodyCoefficients::
    thirdBodyCoefficients(const cu_ushort nsp, const cu_ushort nrc,
                          const cu_ushort nrcThird,
                          const cu_short *__restrict__ indexPtr,
                          const cu_real *__restrict__ coefThirdPtr,
                          const cudaStreams &streams)
    : index_(nrc, indexPtr, streams()), coefThird_(nrcThird, nsp) {
    coefThird_.copyFromHostAsync(coefThirdPtr, streams());
}

inline OpenHurricane::cuChem::thirdBodyCoefficients::
    thirdBodyCoefficients(
        const cu_ushort nsp, const cu_ushort nrc,
        const cuChemInterface::thirdBodyCoeffInterface &thirdBInt)
    : index_(nrc, thirdBInt.thidrBodyIndexPtr_),
      coefThird_(thirdBInt.nrcThird_, nsp, thirdBInt.coefThirdPtr_) {}

inline OpenHurricane::cuChem::thirdBodyCoefficients::
    thirdBodyCoefficients(
        const cu_ushort nsp, const cu_ushort nrc,
        const cuChemInterface::thirdBodyCoeffInterface &thirdBInt,
        const cudaStreams &streams)
    : index_(nrc, thirdBInt.thidrBodyIndexPtr_, streams()),
      coefThird_(thirdBInt.nrcThird_, nsp) {
    coefThird_.copyFromHostAsync(thirdBInt.coefThirdPtr_, streams());
}

inline cu_dual
OpenHurricane::cuChem::thirdBodyCoefficients::thirdBodyCoefficients(
    const thirdBodyCoefficients &tbc)
    : index_(tbc.index_), coefThird_(tbc.coefThird_) {}

inline cu_dual OpenHurricane::cuChem::thirdBodyCoefficients::
    ~thirdBodyCoefficients() noexcept {}

inline cu_host void
OpenHurricane::cuChem::thirdBodyCoefficients::destroy() {
    destroyCuArray(index_);
    destroyCuArray(coefThird_);
}

#endif // CUDA_PARALLEL
