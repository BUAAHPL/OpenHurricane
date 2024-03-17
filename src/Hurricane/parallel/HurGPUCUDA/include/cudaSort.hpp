/*!
 * \file cudaSort.hpp
 * \brief Header of CUDA sort function.
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
#include "CUDAFunctions.hpp"

namespace OpenHurricane {
    namespace cudaSort {
        /**
         * \brief Monolithic Bacther's sort kernel for short arrays.
         */
        template <class keyType, class valueType>
        inline cu_device void
        oddEvenSort(keyType *__restrict__ key, valueType *__restrict__ val, const cu_integer len,
                    const cu_integer tid, const cu_integer dir);

        /**
         * \brief Monolithic Bacther's sort kernel for short arrays.
         */
        template <class keyType>
        inline cu_device void oddEvenSort(keyType *__restrict__ key, const cu_integer len,
                                                 const cu_integer tid, const cu_integer dir);
    } // namespace cudaSort
} // namespace OpenHurricane
#include "cudaSort.inl"

#endif // CUDA_PARALLEL