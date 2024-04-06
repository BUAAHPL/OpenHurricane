#include "cudaReduction.hpp"
/*!
 * \file cudaReduction.inl
 * \brief The In-Line functions of the <i>cudaReduction.hpp</i> file.
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

template <unsigned int dimX>
inline cu_device void OpenHurricane::cudaReduction::warpReduce(volatile cu_real *sPtr,
                                                                  unsigned int tid) {
    if (dimX >= 64) {
        sPtr[tid] += sPtr[tid + 32];
    }
    if (dimX >= 32) {
        if (tid < 16) {
            sPtr[tid] += sPtr[tid + 16];
        }
    }
    if (dimX >= 16) {
        if (tid < 8) {
            sPtr[tid] += sPtr[tid + 8];
        }
    }
    if (dimX >= 8) {
        if (tid < 4) {
            sPtr[tid] += sPtr[tid + 4];
        }
    }
    if (dimX >= 4) {
        if (tid < 2) {
            sPtr[tid] += sPtr[tid + 2];
        }
    }
    if (dimX >= 2) {
        if (tid < 1) {
            sPtr[tid] += sPtr[tid + 1];
        }
    }
}

//template<unsigned int dimX>
//inline cu_device void OpenHurricane::cudaReduction::warpReduce
//(
//	volatile cu_real* sPtr,
//	unsigned int tid,
//	const cu_ushort sSize
//)
//{
//	if (dimX >= 64)
//	{
//		if (tid + 32 < sSize)
//		{
//			sPtr[tid] += sPtr[tid + 32];
//		}
//	}
//	if (dimX >= 32)
//	{
//		if (tid + 16 < sSize)
//		{
//			sPtr[tid] += sPtr[tid + 16];
//		}
//	}
//	if (dimX >= 16)
//	{
//		if (tid + 8 < sSize)
//		{
//			sPtr[tid] += sPtr[tid + 8];
//		}
//	}
//	if (dimX >= 8)
//	{
//		if (tid + 4 < sSize)
//		{
//			sPtr[tid] += sPtr[tid + 4];
//		}
//	}
//	if (dimX >= 4)
//	{
//		if (tid + 2 < sSize)
//		{
//			sPtr[tid] += sPtr[tid + 2];
//		}
//	}
//	if (dimX >= 2)
//	{
//		if (tid + 1 < sSize)
//		{
//			sPtr[tid] += sPtr[tid + 1];
//		}
//	}
//}

template <unsigned int dimX>
inline cu_device void OpenHurricane::cudaReduction::reduce(volatile cu_real *sPtr,
                                                              unsigned int tid) {
    if (dimX >= 512) {
        if (tid < 256) {
            sPtr[tid] += sPtr[tid + 256];
        }
        __syncthreads();
    }
    if (dimX >= 256) {
        if (tid < 128) {
            sPtr[tid] += sPtr[tid + 128];
        }
        __syncthreads();
    }
    if (dimX >= 128) {
        if (tid < 64) {
            sPtr[tid] += sPtr[tid + 64];
        }
        __syncthreads();
    }

    if (tid < 32) {
        warpReduce<dimX>(sPtr, tid);
    }
    __syncthreads();
}

template <unsigned int dimX>
inline cu_device void OpenHurricane::cudaReduction::reduce(volatile cu_real *sPtr,
                                                              unsigned int tid,
                                                              const cu_ushort sSize) {
    auto j = tid + dimX;
    while (tid < sSize && j < sSize) {
        sPtr[tid] += sPtr[j];
        j += dimX;
    }

    __syncthreads();
    reduce<dimX>(sPtr, tid);
}

inline cu_device void OpenHurricane::cudaReduction::warpReduce(volatile cu_real *sPtr,
                                                                  unsigned int tid,
                                                                  const unsigned int dimX) {
    if (dimX >= 64) {
        sPtr[tid] += sPtr[tid + 32];
    }
    if (dimX >= 32) {
        if (tid < 16) {
            sPtr[tid] += sPtr[tid + 16];
        }
    }
    if (dimX >= 16) {
        if (tid < 8) {
            sPtr[tid] += sPtr[tid + 8];
        }
    }
    if (dimX >= 8) {
        if (tid < 4) {
            sPtr[tid] += sPtr[tid + 4];
        }
    }
    if (dimX >= 4) {
        if (tid < 2) {
            sPtr[tid] += sPtr[tid + 2];
        }
    }
    if (dimX >= 2) {
        if (tid < 1) {
            sPtr[tid] += sPtr[tid + 1];
        }
    }
}

inline cu_device void OpenHurricane::cudaReduction::reduce(volatile cu_real *sPtr,
                                                              unsigned int tid,
                                                              const unsigned int dimX) {
    if (dimX >= 512) {
        if (tid < 256) {
            sPtr[tid] += sPtr[tid + 256];
        }
        __syncthreads();
    }
    if (dimX >= 256) {
        if (tid < 128) {
            sPtr[tid] += sPtr[tid + 128];
        }
        __syncthreads();
    }
    if (dimX >= 128) {
        if (tid < 64) {
            sPtr[tid] += sPtr[tid + 64];
        }
        __syncthreads();
    }

    if (tid < 32) {
        warpReduce(sPtr, tid, dimX);
    }
    __syncthreads();
}

cu_device void OpenHurricane::cudaReduction::reduce(volatile cu_real *sPtr, unsigned int tid,
                                                       const cu_ushort sSize,
                                                       const unsigned int dimX) {
    auto j = tid + dimX;
    while (tid < sSize && j < sSize) {
        sPtr[tid] += sPtr[j];
        j += dimX;
    }

    __syncthreads();
    reduce(sPtr, tid, dimX);
}

template <unsigned int dimX>
inline cu_device void OpenHurricane::cudaReduction::warpReduceMax(volatile cu_real *sPtr,
                                                                     unsigned int tid) {
    if (dimX >= 64) {
        sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 32]);
    }
    if (dimX >= 32) {
        if (tid < 16) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 16]);
        }
    }
    if (dimX >= 16) {
        if (tid < 8) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 8]);
        }
    }
    if (dimX >= 8) {
        if (tid < 4) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 4]);
        }
    }
    if (dimX >= 4) {
        if (tid < 2) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 2]);
        }
    }
    if (dimX >= 2) {
        if (tid < 1) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 1]);
        }
    }
}

template <unsigned int dimX>
inline cu_device void OpenHurricane::cudaReduction::reduceMax(volatile cu_real *sPtr,
                                                                 unsigned int tid) {
    if (dimX >= 512) {
        if (tid < 256) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 256]);
        }
        __syncthreads();
    }
    if (dimX >= 256) {
        if (tid < 128) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 128]);
        }
        __syncthreads();
    }
    if (dimX >= 128) {
        if (tid < 64) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 64]);
        }
        __syncthreads();
    }

    if (tid < 32) {
        warpReduceMax<dimX>(sPtr, tid);
    }
    __syncthreads();
}

template <unsigned int dimX>
inline cu_device void OpenHurricane::cudaReduction::reduceMax(volatile cu_real *sPtr,
                                                                 unsigned int tid,
                                                                 const cu_ushort sSize) {
    auto j = tid + dimX;
    while (tid < sSize && j < sSize) {
        sPtr[tid] = cu_max(sPtr[tid], sPtr[j]);
        j += dimX;
    }
    __syncthreads();
    reduceMax<dimX>(sPtr, tid);
}

inline cu_device void OpenHurricane::cudaReduction::warpReduceMax(volatile cu_real *sPtr,
                                                                     unsigned int tid,
                                                                     const unsigned int dimX) {
    if (dimX >= 64) {
        sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 32]);
    }
    if (dimX >= 32) {
        if (tid < 16) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 16]);
        }
    }
    if (dimX >= 16) {
        if (tid < 8) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 8]);
        }
    }
    if (dimX >= 8) {
        if (tid < 4) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 4]);
        }
    }
    if (dimX >= 4) {
        if (tid < 2) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 2]);
        }
    }
    if (dimX >= 2) {
        if (tid < 1) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 1]);
        }
    }
}

inline cu_device void OpenHurricane::cudaReduction::reduceMax(volatile cu_real *sPtr,
                                                                 unsigned int tid,
                                                                 const unsigned int dimX) {
    if (dimX >= 512) {
        if (tid < 256) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 256]);
        }
        __syncthreads();
    }
    if (dimX >= 256) {
        if (tid < 128) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 128]);
        }
        __syncthreads();
    }
    if (dimX >= 128) {
        if (tid < 64) {
            sPtr[tid] = cu_max(sPtr[tid], sPtr[tid + 64]);
        }
        __syncthreads();
    }

    if (tid < 32) {
        warpReduceMax(sPtr, tid, dimX);
    }
    __syncthreads();
}

inline cu_device void OpenHurricane::cudaReduction::reduceMax(volatile cu_real *sPtr,
                                                                 unsigned int tid,
                                                                 const cu_ushort sSize,
                                                                 const unsigned int dimX) {
    auto j = tid + dimX;
    while (tid < sSize && j < sSize) {
        sPtr[tid] = cu_max(sPtr[tid], sPtr[j]);
        j += dimX;
    }
    __syncthreads();
    reduceMax(sPtr, tid, dimX);
}

template <unsigned int dimX>
inline cu_device void OpenHurricane::cudaReduction::warpReduceMin(volatile cu_real *sPtr,
                                                                     unsigned int tid) {
    if (dimX >= 64) {
        sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 32]);
    }
    if (dimX >= 32) {
        if (tid < 16) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 16]);
        }
    }
    if (dimX >= 16) {
        if (tid < 8) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 8]);
        }
    }
    if (dimX >= 8) {
        if (tid < 4) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 4]);
        }
    }
    if (dimX >= 4) {
        if (tid < 2) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 2]);
        }
    }
    if (dimX >= 2) {
        if (tid < 1) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 1]);
        }
    }
}

template <unsigned int dimX>
inline cu_device void OpenHurricane::cudaReduction::reduceMin(volatile cu_real *sPtr,
                                                                 unsigned int tid) {
    if (dimX >= 512) {
        if (tid < 256) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 256]);
        }
        __syncthreads();
    }
    if (dimX >= 256) {
        if (tid < 128) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 128]);
        }
        __syncthreads();
    }
    if (dimX >= 128) {
        if (tid < 64) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 64]);
        }
        __syncthreads();
    }

    if (tid < 32) {
        warpReduceMin<dimX>(sPtr, tid);
    }
    __syncthreads();
}

template <unsigned int dimX>
inline cu_device void OpenHurricane::cudaReduction::reduceMin(volatile cu_real *sPtr,
                                                                 unsigned int tid,
                                                                 const cu_ushort sSize) {
    auto j = tid + dimX;
    while (tid < sSize && j < sSize) {
        sPtr[tid] = cu_min(sPtr[tid], sPtr[j]);
        j += dimX;
    }
    __syncthreads();
    reduceMin<dimX>(sPtr, tid);
}

inline cu_device void OpenHurricane::cudaReduction::warpReduceMin(volatile cu_real *sPtr,
                                                                     unsigned int tid,
                                                                     const unsigned int dimX) {
    if (dimX >= 512) {
        if (tid < 256) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 256]);
        }
        __syncthreads();
    }
    if (dimX >= 256) {
        if (tid < 128) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 128]);
        }
        __syncthreads();
    }
    if (dimX >= 128) {
        if (tid < 64) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 64]);
        }
        __syncthreads();
    }

    if (tid < 32) {
        warpReduceMin(sPtr, tid, dimX);
    }
    __syncthreads();
}

inline cu_device void OpenHurricane::cudaReduction::reduceMin(volatile cu_real *sPtr,
                                                                 unsigned int tid,
                                                                 const unsigned int dimX) {
    if (dimX >= 512) {
        if (tid < 256) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 256]);
        }
        __syncthreads();
    }
    if (dimX >= 256) {
        if (tid < 128) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 128]);
        }
        __syncthreads();
    }
    if (dimX >= 128) {
        if (tid < 64) {
            sPtr[tid] = cu_min(sPtr[tid], sPtr[tid + 64]);
        }
        __syncthreads();
    }

    if (tid < 32) {
        warpReduceMin(sPtr, tid, dimX);
    }
    __syncthreads();
}

inline cu_device void OpenHurricane::cudaReduction::reduceMin(volatile cu_real *sPtr,
                                                                 unsigned int tid,
                                                                 const cu_ushort sSize,
                                                                 const unsigned int dimX) {
    auto j = tid + dimX;
    while (tid < sSize && j < sSize) {
        sPtr[tid] = cu_min(sPtr[tid], sPtr[j]);
        j += dimX;
    }
    __syncthreads();
    reduceMin(sPtr, tid, dimX);
}

#endif // CUDA_PARALLEL
