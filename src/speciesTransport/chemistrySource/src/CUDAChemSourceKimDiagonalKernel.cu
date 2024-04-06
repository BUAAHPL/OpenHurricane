/*!
 * \file CUDAChemSourceKimDiagonalKernel.cu
 * \brief Main subroutines for computing chemistry source with Kim diagonal Jacobian in CUDA platform.
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

#ifdef CUDA_PARALLEL
#include "CUDAFunctions.hpp"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "CUDAChemSourceKimDiagonalKernel.h"
#include <stdio.h>

namespace OpenHurricane {
    template<class dRiType>
    __global__ void calcChemSourceTermsKimDiagonalDevice(cu2DArray<cu_real> yiRhoT,
                                                         const cuChem::reactionTable reactions,
                                                         cu2DArray<cu_real> Ri,
                                                         cu2DArray<dRiType> dRidrhoyi,
                                                         const cu_integer nCells) {
        // Shared memory on GPU
        extern __shared__ cu_real myshare[];

        const auto nsp = reactions.nsp();
        const auto nrc = reactions.nrc();
        cu_real *ciS = myshare;
        cu_real *omegaiS = (cu_real *)&ciS[blockDim.y * nsp];
        cu_real *GdRTS = (cu_real *)&omegaiS[blockDim.y * nsp];

        // The index of this cell
        const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

        // The dimensional density of this cell [kg/m^3]
        cu_real rhod = 0;

        // The dimensional temperature of this cell [K]
        cu_real Td = 0;

        // To get the dimensional density and temperature.
        if (cellId < nCells) {
            rhod = yiRhoT(cellId, nsp);
            Td = yiRhoT(cellId, nsp + 1);
        }

        // To get the pointer of ci (the molar concentrations) of this cell from the allocated shared memory.
        cu_real *ci = (cu_real *)(ciS + threadIdx.y * nsp);
        dRiType *dOmegaidCi = dRidrhoyi.row(cellId);

        // To get the pointer of omegai (the chemical source terms) of this cell from the allocated shared memory.
        cu_real *omegai = (cu_real *)(omegaiS + threadIdx.y * (nsp));

        // To get the pointer of GdRT (the Gibb's free energy) of this cell from the allocated shared memory.
        cu_real *GdRT = (cu_real *)(GdRTS + threadIdx.y * nsp);

        // To get the molar concentrations and the Gibb's free energy for per species.
        auto i = threadIdx.x;
        if (cellId < nCells) {
            while (i < nsp) {
                ci[i] = reactions.species().yiToci(rhod, yiRhoT(cellId, i), i);
                omegai[i] = 0;
                GdRT[i] = reactions.species().nGdRT(i, Td);
                if (i < nsp - 1) {
                    dOmegaidCi[i] = 0;
                }
                i += blockDim.x;
            }
        }
        __syncthreads();

        // To get the molar chemical source terms from each reaction.
        i = threadIdx.x;
        if (cellId < nCells) {
            while (i < nrc) {
                reactions.omegaCoupled3<dRiType>(i, Td, ci, GdRT, omegai, dOmegaidCi);
                i += blockDim.x;
            }
        }
        __syncthreads();

        // To get the molar chemical source terms from each reaction.
        i = threadIdx.x;
        if (cellId < nCells) {
            while (i < nsp) {
                reactions.Ri(i, omegai, Ri.row(cellId));

                i += blockDim.x;
            }
        }
    }

    void getBlockAndGridSize1(const cu_integer nsp, const cu_integer nCells, dim3 &blockSize,
                              dim3 &gridSize) {
        unsigned int blocksize_x, blocksize_y;
        unsigned int gridsize_x;
        if (nsp < 8) {
            blocksize_x = 4;
            blocksize_y = 32;
            //blocksize_y = 64;
            gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
        } else if (nsp < 16) {
            blocksize_x = 8;
            blocksize_y = 16;
            //blocksize_y = 32;
            gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
        } else if (nsp < 32) {
            blocksize_x = 16;
            blocksize_y = 8;
            //blocksize_y = 16;
            gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
        } else if (nsp < 64) {
            blocksize_x = 32;
            blocksize_y = 4;
            //blocksize_y = 8;
            gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
        } else if (nsp < 128) {
            blocksize_x = 64;
            blocksize_y = 2;
            //blocksize_y = 4;
            gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
        } else {
            blocksize_x = 128;
            blocksize_y = 1;
            //blocksize_y = 2;
            gridsize_x = (nCells + blocksize_y - 1) / blocksize_y;
        }

        dim3 blockSizeTmp(blocksize_x, blocksize_y, 1);
        dim3 gridSizeTmp(gridsize_x, 1, 1);

        blockSize = blockSizeTmp;
        gridSize = gridSizeTmp;
    }

    size_t getSharedMemSize1(const cu_integer nsp, const unsigned int blocksize_y) {
        return blocksize_y * (nsp + nsp + nsp) * sizeof(cu_real);
    }

    void onlyChemistrySource::calcChemSourceTermsKimDiagonal(
        const cu_real *__restrict__ hostYiRhoT, const cuChem::reactionTable &reactions,
        const cu_integer nsp, const cu_integer nrc, const cu_integer nCells,
        cu_real *__restrict__ Ri, cu_real *__restrict__ dRidrhoyi) {
        cu2DArray<cu_real> dYiRhoT(nsp + 2, nCells, hostYiRhoT);
        cu2DArray<cu_real> dRi(nsp, nCells);
        cu2DArray<cu_real> device_dRidrhoyi(nsp - 1, nCells);

        dim3 blockSize;
        dim3 gridSize;

        getBlockAndGridSize1(nsp, nCells, blockSize, gridSize);

        const auto shareMemSize = getSharedMemSize1(nsp, blockSize.y);
        calcChemSourceTermsKimDiagonalDevice<<<gridSize, blockSize, shareMemSize, 0>>>(
            dYiRhoT, reactions, dRi, device_dRidrhoyi, nCells);

        checkCUDAError(cudaGetLastError());
        dYiRhoT.clear();
        dRi.copyToHost(Ri);
        device_dRidrhoyi.copyToHost(dRidrhoyi);
        dRi.clear();
        device_dRidrhoyi.clear();
    }
    void onlyChemistrySource::calcChemSourceTermsKimDiagonal(
        const cu_real *__restrict__ hostYiRhoT, const cuChem::reactionTable &reactions,
        const cu_integer nsp, const cu_integer nrc, const cu_integer nCells,
        const nThreadsAndBlocks &nTB, cu_real *__restrict__ Ri,
        cu_real *__restrict__ dRidrhoyi) {
        cu2DArray<cu_real> dYiRhoT(nsp + 2, nCells, hostYiRhoT);
        cu2DArray<cu_real> dRi(nsp, nCells);
        cu2DArray<cu_real> device_dRidrhoyi(nsp - 1, nCells);
        calcChemSourceTermsKimDiagonalDevice<<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(),
                                               0>>>(dYiRhoT, reactions, dRi, device_dRidrhoyi,
                                                    nCells);

        checkCUDAError(cudaGetLastError());
        dYiRhoT.clear();
        dRi.copyToHost(Ri);
        device_dRidrhoyi.copyToHost(dRidrhoyi);
        dRi.clear();
        device_dRidrhoyi.clear();
    }

    void onlyChemistrySource::calcChemSourceTermsKimDiagonalAsync(
        const cu_real *__restrict__ hostYiRhoT, cu2DArray<cu_real> &dYiRhoT,
        const cuChem::reactionTable &reactions, const cu_integer nsp, const cu_integer nrc,
        const cu_integer nCells, const nThreadsAndBlocks &nTB, cu_real *__restrict__ Ri,
        cu_real *__restrict__ dRidrhoyi, cu2DArray<cu_real> &dRi,
        cu2DArray<cu_real> &ddRidrhoyi, const cudaStreams &streams) {
        dYiRhoT.copyFromHostAsync(hostYiRhoT, streams());
        calcChemSourceTermsKimDiagonalDevice<<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(),
                                               streams()>>>(dYiRhoT, reactions, dRi, ddRidrhoyi,
                                                            nCells);

        dRi.copyToHostAsync(Ri, streams());
        ddRidrhoyi.copyToHostAsync(dRidrhoyi, streams());
    }

    void onlyChemistrySource::calcChemSourceTermsKimDiagonalAsyncHybrid(
        const cu_real *__restrict__ hostYiRhoT, cu2DArray<cu_real> &dYiRhoT,
        const cuChem::reactionTable &reactions, const cu_integer nsp, const cu_integer nrc,
        const cu_integer nCells, const nThreadsAndBlocks &nTB, cu_real *__restrict__ Ri,
        cu_float *__restrict__ dRidrhoyi, cu2DArray<cu_real> &dRi,
        cu2DArray<cu_float> &ddRidrhoyi, const cudaStreams &streams) {
        dYiRhoT.copyFromHostAsync(hostYiRhoT, streams());
        calcChemSourceTermsKimDiagonalDevice<cu_float>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(),
                                               streams()>>>(dYiRhoT, reactions, dRi, ddRidrhoyi,
                                                            nCells);

        dRi.copyToHostAsync(Ri, streams());
        ddRidrhoyi.copyToHostAsync(dRidrhoyi, streams());
    }

} // namespace OpenHurricane

#endif // CUDA_PARALLEL