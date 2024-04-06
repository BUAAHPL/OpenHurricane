
/*!
 * \file cuAttributes.cu
 * \brief Subroutines of the <i>cuAttributes.hpp</i> file.
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

#include "cuAttributes.hpp"
#include "cuPreset.hpp"
#ifdef CUDA_PARALLEL

bool OpenHurricane::cuAttributes::deviceSupportsMemoryPools_ = false;
bool OpenHurricane::cuAttributes::deviceSupportsConcurrentManagedAccess_ = false;
int OpenHurricane::cuAttributes::cudaDriverVersion_ = 0;

void OpenHurricane::cuAttributes::init() {
    checkCUDAError(cudaDriverGetVersion(&cudaDriverVersion_));
    if (cudaDriverVersion_ >= 11020) {
        deviceSupportsMemoryPools_ = true;
        deviceSupportsConcurrentManagedAccess_ = true;
        int deviceCount;

        checkCUDAError(cudaGetDeviceCount(&deviceCount));
        for (int ig = 0; ig < deviceCount; ++ig) {
            int val;
            checkCUDAError(cudaDeviceGetAttribute(&val, cudaDevAttrMemoryPoolsSupported, ig));
            if (val == 0) {
                deviceSupportsMemoryPools_ = false;
                break;
            }

            checkCUDAError(cudaDeviceGetAttribute(&val, cudaDevAttrConcurrentManagedAccess, ig));
            if (val == 0) {
                deviceSupportsConcurrentManagedAccess_ = false;
                break;
            }
        }
    }
}

#else // CUDA_PARALLEL

#endif // CUDA_PARALLEL
