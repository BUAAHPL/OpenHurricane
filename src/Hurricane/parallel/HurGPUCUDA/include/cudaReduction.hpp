/*!
 * \file cudaReduction.hpp
 * \brief Header of CUDA reduction.
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
#include "CUDAFunctions.hpp"

namespace OpenHurricane {
    namespace cudaReduction {
        // +

        template <unsigned int dimX>
        inline cu_device void warpReduce(volatile cu_real *sPtr, unsigned int tid);

        template <unsigned int dimX>
        inline cu_device void reduce(volatile cu_real *sPtr, unsigned int tid);

        template <unsigned int dimX>
        inline cu_device void reduce(volatile cu_real *sPtr, unsigned int tid,
                                            const cu_ushort sSize);

        inline cu_device void warpReduce(volatile cu_real *sPtr, unsigned int tid,
                                                const unsigned int dimX);

        inline cu_device void reduce(volatile cu_real *sPtr, unsigned int tid,
                                            const unsigned int dimX);

        inline cu_device void reduce(volatile cu_real *sPtr, unsigned int tid,
                                            const cu_ushort sSize, const unsigned int dimX);

        // max()

        template <unsigned int dimX>
        inline cu_device void warpReduceMax(volatile cu_real *sPtr, unsigned int tid);

        template <unsigned int dimX>
        inline cu_device void reduceMax(volatile cu_real *sPtr, unsigned int tid);

        template <unsigned int dimX>
        inline cu_device void reduceMax(volatile cu_real *sPtr, unsigned int tid,
                                               const cu_ushort sSize);

        inline cu_device void warpReduceMax(volatile cu_real *sPtr, unsigned int tid,
                                                   const unsigned int dimX);

        inline cu_device void reduceMax(volatile cu_real *sPtr, unsigned int tid,
                                               const unsigned int dimX);

        inline cu_device void reduceMax(volatile cu_real *sPtr, unsigned int tid,
                                               const cu_ushort sSize, const unsigned int dimX);

        // min()

        template <unsigned int dimX>
        inline cu_device void warpReduceMin(volatile cu_real *sPtr, unsigned int tid);

        template <unsigned int dimX>
        inline cu_device void reduceMin(volatile cu_real *sPtr, unsigned int tid);

        template <unsigned int dimX>
        inline cu_device void reduceMin(volatile cu_real *sPtr, unsigned int tid,
                                               const cu_ushort sSize);

        inline cu_device void warpReduceMin(volatile cu_real *sPtr, unsigned int tid,
                                                   const unsigned int dimX);

        inline cu_device void reduceMin(volatile cu_real *sPtr, unsigned int tid,
                                               const unsigned int dimX);

        inline cu_device void reduceMin(volatile cu_real *sPtr, unsigned int tid,
                                               const cu_ushort sSize, const unsigned int dimX);
    } // namespace cudaReduction
} // namespace OpenHurricane
#include "cudaReduction.inl"

#endif // CUDA_PARALLEL