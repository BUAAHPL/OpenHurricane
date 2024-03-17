/*!
 * \file thermoListCUDA.cu
 * \brief The subroutines and functions of thermo in CUDA platform.
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

#include "thermoListCUDA.hpp"
#include <cmath>
#ifdef CUDA_PARALLEL

#include "cudaReduction.hpp"

namespace OpenHurricane {
    namespace CUDAThermo {
        __global__ void gethai(cu2DArray<cu_real> yiTP, const cuChem::speciesTable spcs,
                               cu2DArray<cu_real> hai, const cu_integer nsp,
                               const cu_integer nCells) {
            // The index of this cell
            const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

            // The dimensional temperature of this cell [K]
            cu_real Td = 0;

            // To get the dimensional density and temperature.
            if (cellId < nCells) {
                Td = yiTP(cellId, nsp);
            }

            auto i = threadIdx.x;
            if (cellId < nCells) {
                while (i < nsp) {
                    hai(cellId, i) = spcs.ha0(Td, i);
                    i += blockDim.x;
                }
            }
            __syncthreads();
        }

        void calchai(const cu_real *__restrict__ hostYiTP, const cuChem::speciesTable &spcs,
                     const nThreadsAndBlocks &nTB, const cu_integer nsp,
                     const cu_integer nCells, cu_real *__restrict__ host_hai) {
            cu2DArray<cu_real> dYiTP(nsp + 2, nCells, hostYiTP);
            cu2DArray<cu_real> d_hai(nsp, nCells);

            gethai<<<nTB.gridSize(), nTB.blockSize(), 0>>>(dYiTP, spcs, d_hai, nsp, nCells);
            d_hai.copyToHost(host_hai);
        }

        void calchaiAsync(const cu_real *__restrict__ hostYiTP, cu2DArray<cu_real> &dYiTP,
                          const cuChem::speciesTable &spcs, const nThreadsAndBlocks &nTB,
                          const cu_integer nsp, const cu_integer nCells,
                          cu_real *__restrict__ host_hai, cu2DArray<cu_real> &d_hai,
                          const cudaStreams &streams) {
            dYiTP.copyFromHostAsync(hostYiTP, streams());
            gethai<<<nTB.gridSize(), nTB.blockSize(), 0, streams()>>>(dYiTP, spcs, d_hai, nsp,
                                                                      nCells);
            d_hai.copyToHostAsync(host_hai, streams());
        }

        template <unsigned int dimX>
        __global__ void gethaicp(cu2DArray<cu_real> yiTP, const cuChem::speciesTable spcs,
                                 cu2DArray<cu_real> haicp, const cu_integer nsp,
                                 const cu_integer nCells) {
            extern __shared__ cu_real myshare[];

            cu_real *cpis = myshare;

            // The index of this cell
            const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

            // The dimensional temperature of this cell [K]
            cu_real Td = 0;

            // To get the dimensional density and temperature.
            if (cellId < nCells) {
                Td = yiTP(cellId, nsp);
            }

            cu_real *cpi = (cu_real *)(cpis + threadIdx.y * nsp);

            auto i = threadIdx.x;
            if (cellId < nCells) {
                while (i < nsp) {
                    cpi[i] = spcs.cp0(Td, i) * yiTP(cellId, i);
                    haicp(cellId, i) = spcs.ha0(Td, i);
                    i += blockDim.x;
                }
            }
            __syncthreads();

            i = threadIdx.x;
            cudaReduction::reduce<dimX>(cpi, i, nsp);

            i = threadIdx.x;
            if (cellId < nCells) {
                if (i == 0) {
                    haicp(cellId, nsp) = cpi[0];
                }
            }
        }

        void calchaicp(const cu_real *__restrict__ hostYiTP, const cuChem::speciesTable &spcs,
                       const nThreadsAndBlocks &nTB, const cu_integer nsp,
                       const cu_integer nCells, cu_real *__restrict__ host_haicp) {
            cu2DArray<cu_real> dYiTP(nsp + 2, nCells, hostYiTP);
            cu2DArray<cu_real> d_haicp(nsp + 1, nCells);

            if (nTB.blockSize().x == 4) {
                gethaicp<4><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else if (nTB.blockSize().x == 8) {
                gethaicp<8><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else if (nTB.blockSize().x == 16) {
                gethaicp<16><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else if (nTB.blockSize().x == 32) {
                gethaicp<32><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else if (nTB.blockSize().x == 64) {
                gethaicp<64><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else if (nTB.blockSize().x == 128) {
                gethaicp<128><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else {
                fprintf(stderr,
                        "OpenHurricane Error Of Calling CUDA function: Unsupported block "
                        "size in x direction: %d at %s:%d of \"%s\" \n",
                        nTB.blockSize().x, HUR_FILE, HUR_LINE, HUR_FUNCTION);
                exit(EXIT_FAILURE);
            }
            d_haicp.copyToHost(host_haicp);
        }

        void calchaicpAsync(const cu_real *__restrict__ hostYiTP, cu2DArray<cu_real> &dYiTP,
                            const cuChem::speciesTable &spcs, const nThreadsAndBlocks &nTB,
                            const cu_integer nsp, const cu_integer nCells,
                            cu_real *__restrict__ host_haicp, cu2DArray<cu_real> &d_haicp,
                            const cudaStreams &streams) {
            dYiTP.copyFromHostAsync(hostYiTP, streams());
            if (nTB.blockSize().x == 4) {
                gethaicp<4><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), streams()>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else if (nTB.blockSize().x == 8) {
                gethaicp<8><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), streams()>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else if (nTB.blockSize().x == 16) {
                gethaicp<16><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), streams()>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else if (nTB.blockSize().x == 32) {
                gethaicp<32><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), streams()>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else if (nTB.blockSize().x == 64) {
                gethaicp<64><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), streams()>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else if (nTB.blockSize().x == 128) {
                gethaicp<128><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), streams()>>>(
                    dYiTP, spcs, d_haicp, nsp, nCells);
            } else {
                fprintf(stderr,
                        "OpenHurricane Error Of Calling CUDA function: Unsupported block "
                        "size in x direction: %d at %s:%d of \"%s\" \n",
                        nTB.blockSize().x, HUR_FILE, HUR_LINE, HUR_FUNCTION);
                exit(EXIT_FAILURE);
            }
            d_haicp.copyToHostAsync(host_haicp, streams());
        }

    } // namespace CUDAThermo
} // namespace OpenHurricane

#endif // CUDA_PARALLEL