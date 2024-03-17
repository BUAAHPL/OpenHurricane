/*!
 * \file GPUInfo.hpp
 * \brief Headers of the GPU information.
 *        The subroutines and functions are in the <i>GPUInfo.cpp</i> file.
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
#include "logFile.hpp"

#ifdef CUDA_PARALLEL
#if CUDART_VERSION < 5000
// CUDA-C includes
#include <cuda.h>

// This function wraps the CUDA Driver API into a template function
template <class T>
inline void getCudaAttribute(T *attribute, CUdevice_attribute device_attribute, int device) {
    CUresult error = cuDeviceGetAttribute(attribute, device_attribute, device);

    if (CUDA_SUCCESS != error) {
        fprintf(stderr,
                "cuSafeCallNoSync() Driver API error = %04d from file <%s>, "
                "line %i.\n",
                error, __FILE__, __LINE__);

        exit(EXIT_FAILURE);
    }
}
#endif // CUDART_VERSION < 5000
#endif // CUDA_PARALLEL

namespace OpenHurricane {
    class GPUInfo {
    public:
        // Static functions

        /*!\brief To get the information of GPU devices.
                Return 0 if not have GPU device that support CUDA.
                Return the count of the device if succeed.*/
        static short int getGPUInfo();

        static int getGPUNum();
    };

} // namespace OpenHurricane
