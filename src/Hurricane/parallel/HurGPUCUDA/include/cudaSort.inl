#include "cudaSort.hpp"
/*!
 * \file cudaSort.inl
 * \brief The In-Line functions of the <i>cudaSort.hpp</i> file.
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

template <class keyType, class valueType>
inline cu_device void
OpenHurricane::cudaSort::oddEvenSort(keyType *__restrict__ key, valueType *__restrict__ val,
                                 const cu_integer len, const cu_integer tid,
                                 const cu_integer dir) {
    unsigned short isOdd = 1;
    auto count = len;
    do {
        __syncthreads();
        isOdd = (++isOdd) % 2;
        const cu_integer id1 = isOdd + 2 * tid;
        if (id1 >= len) {
            continue;
        }
        const cu_integer id2 = isOdd + 2 * tid + 1;
        if (id2 >= len) {
            continue;
        }
        if ((key[id1] > key[id2]) == dir) {
            cuSwap(key[id1], key[id2]);
            cuSwap(val[id1], val[id2]);
        }
    } while (count--);
}

template <class keyType>
inline cu_device void
OpenHurricane::cudaSort::oddEvenSort(keyType *__restrict__ key, const cu_integer len,
                                 const cu_integer tid, const cu_integer dir) {
    unsigned short isOdd = 1;
    auto count = len;
    do {
        __syncthreads();
        isOdd = (++isOdd) % 2;
        const cu_integer id1 = isOdd + 2 * tid;
        if (id1 >= len) {
            continue;
        }
        const cu_integer id2 = isOdd + 2 * tid + 1;
        if (id2 >= len) {
            continue;
        }
        if ((key[id1] > key[id2]) == dir) {
            cuSwap(key[id1], key[id2]);
        }
    } while (count--);
}

#endif // CUDA_PARALLEL
