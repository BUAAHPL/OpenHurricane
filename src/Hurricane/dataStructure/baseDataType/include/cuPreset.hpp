/*!
 * \file cuPreset.hpp
 * \brief Headers of pre-set code for CUDA.
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

#include "preset.hpp"

#ifdef CUDA_PARALLEL
#include <cuda_runtime.h>
#include <device_atomic_functions.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <string>

#ifdef cu_device
#undef cu_device
#endif // cu_device

#ifdef cu_host
#undef cu_host
#endif // cu_host

// CUDA device (GPU) function flag
#define cu_device __device__
#define cu_global __global__
// CUDA host (CPU) function flag
#define cu_host __host__
// CUDA host (CPU) and device (GPU) function flag
#define cu_dual cu_device cu_host

#ifdef HURRICANE_32_INT
using cu_integer = int32_t;
using cu_uinteger = uint32_t;
#elif defined HURRICANE_64_INT
using cu_integer = int64_t;
using cu_uinteger = uint64_t;
#else
using cu_integer = int32_t;
using cu_uinteger = uint32_t;
#endif // HURRICANE_32_INT

#if defined(HURRICANE_SP)
using cu_real = float;
using cu_float = float;
#define cu_tiny 1.0e-6f
#define cu_veryTiny 1.0e-37f
#define cu_rootLarge 1.0e+18f
#elif defined(HURRICANE_DP)
using cu_real = double;
using cu_float = float;
#define cu_tiny 1.0e-15
#define cu_veryTiny 1.0e-300
#define cu_rootLarge 1.0e+150
#endif

using cu_short = short;
using cu_ushort = unsigned short;

#ifndef cu_max
#define cu_max(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef cu_min
#define cu_min(a, b) (((a) < (b)) ? (a) : (b))
#endif

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
//nothing to be done
#else
//atomicAdd() for double-precision floating-point numbers is not available
//on devices with compute capability lower than 6.0 but it can be implemented as follows:
static __inline__ cu_device double atomicAdd(double *address, double val) {
    unsigned long long int *address_as_ull = (unsigned long long int *)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

// Beginning of GPU Architecture definitions
hur_nodiscard inline int _ConvertSMVer2Cores(int major, int minor) {
    // Defines for GPU Architecture types (using the SM version to determine
    // the # of cores per SM
    typedef struct {
        int SM; // 0xMm (hexidecimal notation), M = SM Major version,
        // and m = SM minor version
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] = {
        {0x30, 192}, {0x32, 192}, {0x35, 192}, {0x37, 192}, {0x50, 128}, {0x52, 128},
        {0x53, 128}, {0x60, 64},  {0x61, 128}, {0x62, 128}, {0x70, 64},  {0x72, 64},
        {0x75, 64},  {0x80, 64},  {0x86, 128}, {0x87, 128}, {-1, -1}};

    int index = 0;

    while (nGpuArchCoresPerSM[index].SM != -1) {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
            return nGpuArchCoresPerSM[index].Cores;
        }

        index++;
    }

    // If we don't find the values, we default use the previous one
    // to run properly
    printf("MapSMtoCores for SM %d.%d is undefined."
           "  Default to use %d Cores/SM\n",
           major, minor, nGpuArchCoresPerSM[index - 1].Cores);
    return nGpuArchCoresPerSM[index - 1].Cores;
}

/**
 * \brief If cuda function called succeeded.
 * \param[in] error - Error code
 * \param[out] The description string for the error code
 * \return Return if succeeded.
 */
inline bool isCUDASuccess(const cudaError_t error, std::string &errMsg) {
    errMsg = cudaGetErrorString(error);
    return error == cudaSuccess;
}

template <typename T>
inline void check(T result, char const *const func, const char *const file, int const line) {
    if (result) {
#if defined(HUR_FULL_LOGGER)
        fprintf(stderr, "    Error: CUDA error at %s:%d code=%d(%s) in function: \"%s\" \n", file, line,
                static_cast<unsigned int>(result), cudaGetErrorName(result), func);
#elif defined(HUR_LESS_LOGGER)
        fprintf(stderr, "    Error: CUDA error code=%d(%s) in function: \"%s\" \n", 
                static_cast<unsigned int>(result), cudaGetErrorName(result), func);
#else
        fprintf(stderr, "    Error: CUDA error code=%d(%s)\n", static_cast<unsigned int>(result),
                cudaGetErrorName(result));
#endif // DEBUG
        exit(EXIT_FAILURE);
    }
}

#ifndef checkCUDAError
#define checkCUDAError(val) check((val), #val, HUR_FILE, HUR_LINE)
#endif // !checkCUDAError

#else // CUDA_PARALLEL
// OpenHurricane device (GPU) function flag
#define cu_device
// OpenHurricane host (CPU) function flag
#define cu_host
// OpenHurricane host (CPU) and device (GPU) function flag
#define cu_dual cu_device cu_host
#endif // CUDA_PARALLEL
