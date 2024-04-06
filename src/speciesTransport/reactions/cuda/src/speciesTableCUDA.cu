/*!
 * \file speciesTableCUDA.cu
 * \brief The subroutines and functions of species Table in CUDA platform.
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

#include "speciesTableCUDA.hpp"
#ifdef CUDA_PARALLEL

#include "cudaReduction.hpp"
#include <cmath>

namespace OpenHurricane {
namespace cuChem {
/**
 * \brief Get mixture cp on GPU.
 * \tparam dimX - The The block dimension size in x-direction
 * \note dimX must be satisfied as nsp < 2 x dimX
 */
template <unsigned int dimX>
__global__ void getCpm(cu2DArray<cu_real> yiRhoT,
                       const speciesTable species, cu1DArray<cu_real> cpm,
                       const cu_integer nsp, const cu_integer nCells) {
    extern __shared__ cu_real myshare[];

    // The index of this cell
    const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

    // The dimensional temperature of this cell [K]
    cu_real Td = 0;

    // To get the dimensional density and temperature.
    if (cellId < nCells) {
        Td = yiRhoT(cellId, nsp + 1);
    }

    // To get the pointer of cpi of this cell from the allocated shared memory.
    cu_real *cpi = (cu_real *)(myshare + threadIdx.y * nsp);

    auto i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            cpi[i] = species.cp0(Td, i) * yiRhoT(cellId, i);
            i += blockDim.x;
        }
    }
    __syncthreads();

    i = threadIdx.x;
    cudaReduction::reduce<dimX>(cpi, i, nsp);

    i = threadIdx.x;
    if (cellId < nCells) {
        if (i == 0) {
            cpm(cellId) = cpi[i];
        }
    }
}

__global__ void getHai(cu2DArray<cu_real> yiRhoT,
                       const speciesTable species, cu2DArray<cu_real> hai,
                       const cu_integer nsp, const cu_integer nCells) {
    // The index of this cell
    const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

    // The dimensional temperature of this cell [K]
    cu_real Td = 0;

    // To get the dimensional density and temperature.
    if (cellId < nCells) {
        Td = yiRhoT(cellId, nsp + 1);
    }

    auto i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            hai(cellId, i) = species.ha0(Td, i) * yiRhoT(cellId, i);
            i += blockDim.x;
        }
    }
}

/**
 * \brief Get molar fractions.
 * \tparam dimX - The The block dimension size in x-direction
 * \note dimX must be satisfied as nsp < 2 x dimX
 */
template <unsigned int dimX>
__global__ void getXi(cu2DArray<cu_real> yiRhoT, const speciesTable species,
                      const cu_integer nsp, const cu_integer nCells) {
    extern __shared__ cu_real myshare[];

    // The index of this cell
    const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

    // To get the pointer of cpi of this cell from the allocated shared memory.
    cu_real *xi = (cu_real *)(myshare + threadIdx.y * nsp);

    auto i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            xi[i] = yiRhoT(cellId, i) / species.Wi(i);
            i += blockDim.x;
        }
    }
    __syncthreads();

    i = threadIdx.x;
    cudaReduction::reduce<dimX>(xi, i, nsp);

    const auto Wm = cu_real(1) / xi[0];

    i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            yiRhoT(cellId, i) = yiRhoT(cellId, i) * Wm / species.Wi(i);
            i += blockDim.x;
        }
    }
    __syncthreads();
}
} // namespace CUDAReactions
} // namespace OpenHurricane

#endif // CUDA_PARALLEL