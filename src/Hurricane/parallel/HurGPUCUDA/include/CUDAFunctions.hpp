/*!
 * \file CUDAFunctions.hpp
 * \brief Headers of the main subroutines of the CUDA Functions.
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
#include "cuPreset.hpp"
#include "cuArrayBase.hpp"

#ifdef CUDA_PARALLEL

namespace OpenHurricane {
    
    template <class Type> inline cu_dual void cuSwap(Type &v1, Type &v2) noexcept {
        Type tmp = v1;
        v1 = v2;
        v2 = tmp;
    }
} // namespace OpenHurricane

inline cu_dual cu_real cuSqr(const cu_real s) noexcept {
    return s * s;
}

#endif // CUDA_PARALLEL
