#include "transportListCUDA.hpp"
/*!
 * \file transportListCUDA.cu
 * \brief The subroutines and functions of transport in CUDA platform.
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

#include "transportListCUDA.hpp"
#include <cmath>
#ifdef CUDA_PARALLEL
#include "cudaReduction.hpp"

cu_host
OpenHurricane::CUDATransport::transportListCUDA::transportListCUDA(
    const cu_ushort nsp, const cu_real *__restrict__ ekb,
    const cu_real *__restrict__ sigma,
    const cuChem::speciesTable species)
    : transport_(nsp, ekb, sigma, species) {}

OpenHurricane::CUDATransport::transportListCUDA::transportListCUDA(
    const cu_ushort nsp, const cu_real *__restrict__ ekb,
    const cu_real *__restrict__ sigma,
    const cuChem::speciesTable species, const cudaStreams &streams)
    : transport_(nsp, ekb, sigma, species, streams) {}

cu_device void
OpenHurricane::CUDATransport::transportListCUDA::WilkeCoeff(
    const cu_real *__restrict__ mui, const cu_real *__restrict__ xi,
    const cu_ushort i, cu_real *__restrict__ wci) const {
    cu_real xp = 0;
    const cu_real muii = mui[i];
    const auto Wii = species().Wi(i);
    for (cu_ushort j = 0; j < nsp(); ++j) {
        xp += xi[j] * PhiMuij(muii, mui[j], Wii, species().Wi(j));
    }

    wci[i] = xi[i] / xp;
}

namespace OpenHurricane {
namespace CUDATransport {

/**
 * \brief Get transport properties.
 * \tparam dimX - The The block dimension size in x-direction
 * \note dimX must be satisfied as nsp < 2 x dimX
 */
template <unsigned int dimX>
__global__ void getDim(cu2DArray<cu_real> xiTP,
                       const transportListCUDA tran,
                       cu2DArray<cu_real> Dimm, const cu_integer nsp,
                       const cu_integer nCells) {
    // The index of this cell
    const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

    // The dimensional temperature of this cell [K]
    cu_real Td = 0;

    // The dimensional pressure of this cell [Pa]
    cu_real pd = 0;

    // To get the dimensional density and temperature.
    if (cellId < nCells) {
        Td = xiTP(cellId, nsp);
        pd = xiTP(cellId, nsp + 1);
    }

    auto i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            Dimm(cellId, i) = tran.transport().Dim(Td, pd, xiTP.row(cellId), i);

            i += blockDim.x;
        }
    }
}

/**
 * \brief Get transport properties.
 * \tparam dimX - The The block dimension size in x-direction
 * \note dimX must be satisfied as nsp < 2 x dimX.
 *       And yiT would be changed into xiT
 */
template <unsigned int dimX>
__global__ void
getMuKappaDim(cu2DArray<cu_real> yiTP, const transportListCUDA tran,
              cu2DArray<cu_real> DimmMuKappa, const cu_integer nsp,
              const cu_integer nCells) {
    extern __shared__ cu_real myshare[];

    cu_real *muS = myshare;
    cu_real *kappaS = (cu_real *)&muS[blockDim.y * nsp];

    // The index of this cell
    const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;

    // The dimensional temperature of this cell [K]
    cu_real Td = 0;
    // The dimensional pressure of this cell [Pa]
    cu_real pd = 0;

    // To get the dimensional density and temperature.
    if (cellId < nCells) {
        Td = yiTP(cellId, nsp);
        pd = yiTP(cellId, nsp + 1);
    }

    cu_real *mui = (cu_real *)(muS + threadIdx.y * nsp);
    cu_real *kappai = (cu_real *)(kappaS + threadIdx.y * (nsp));

    //************************************************************************
    // Change yi to xi
    auto i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            tran.species().yidWi(i, yiTP.row(cellId), mui);
            i += blockDim.x;
        }
    }
    __syncthreads();

    i = threadIdx.x;
    cudaReduction::reduce<dimX>(mui, i, nsp);
    const auto wm = cu_real(1) / mui[0];
    i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            tran.species().yi2xi(i, yiTP.row(cellId), wm, mui);
            yiTP.row(cellId)[i] = mui[i];
            i += blockDim.x;
        }
    }
    __syncthreads();
    /*i = threadIdx.x;
    cudaReduction::reduce<dimX>
            (
                    mui,
                    i,
                    nsp
                    );
    i = threadIdx.x;
    if (cellId < nCells)
    {
            while (i < nsp)
            {
                    yiTP.row(cellId)[i] / cu_max(mui[0], cu_veryTiny);
                    i += blockDim.x;
            }
    }
    __syncthreads();*/
    // End of changing yi to xi
    //======================================================================

    i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            mui[i] = tran.transport().mu(Td, i);
            kappai[i] = tran.transport().kappa(
                Td, mui[i], tran.transport().species().cp0(Td, i), i);
            i += blockDim.x;
        }
    }
    __syncthreads();

    i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            tran.WilkeCoeff(mui, yiTP.row(cellId), i, DimmMuKappa.row(cellId));

            i += blockDim.x;
        }
    }
    __syncthreads();

    i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            tran.muKappa(mui, kappai, yiTP.row(cellId), DimmMuKappa.row(cellId),
                         i);

            i += blockDim.x;
        }
    }
    __syncthreads();

    i = threadIdx.x;
    cudaReduction::reduce<dimX>(mui, i, nsp);

    i = threadIdx.x;
    cudaReduction::reduce<dimX>(kappai, i, nsp);

    i = threadIdx.x;
    if (cellId < nCells) {
        if (i == 0) {
            DimmMuKappa(cellId, nsp) = mui[0];
            DimmMuKappa(cellId, nsp + 1) = kappai[0];
        }
    }

    //***************************************************************************
    // Compute mass diffusion coefficients
    i = threadIdx.x;
    if (cellId < nCells) {
        while (i < nsp) {
            DimmMuKappa(cellId, i) =
                tran.transport().Dim(Td, pd, yiTP.row(cellId), i);
            i += blockDim.x;
        }
    }
    //__syncthreads();
}

void calcTransportProperties(const cu_real *__restrict__ hostYiTP,
                             const transportListCUDA &tran,
                             const nThreadsAndBlocks &nTB,
                             const cu_integer nsp, const cu_integer nCells,
                             cu_real *__restrict__ hostDimMuKappa) {
    cu2DArray<cu_real> dYiTP(nsp + 2, nCells, hostYiTP);
    cu2DArray<cu_real> dDimMuKappa(nsp + 2, nCells);

    if (nTB.blockSize().x == 4) {
        getMuKappaDim<4>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else if (nTB.blockSize().x == 8) {
        getMuKappaDim<8>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else if (nTB.blockSize().x == 16) {
        getMuKappaDim<16>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else if (nTB.blockSize().x == 32) {
        getMuKappaDim<32>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else if (nTB.blockSize().x == 64) {
        getMuKappaDim<64>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else if (nTB.blockSize().x == 128) {
        getMuKappaDim<128>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), 0>>>(
                dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else {
        fprintf(stderr,
                "OpenHurricane Error Of Calling CUDA function: Unsupported block "
                "size in x direction: %d at %s:%d of \"%s\" \n",
                nTB.blockSize().x, HUR_FILE, HUR_LINE, HUR_FUNCTION);
        exit(EXIT_FAILURE);
    }
    dDimMuKappa.copyToHost(hostDimMuKappa);
}

void calcTransportPropertiesAsync(
    const cu_real *__restrict__ hostYiTP, cu2DArray<cu_real> &dYiTP,
    const transportListCUDA &tran, const nThreadsAndBlocks &nTB,
    const cu_integer nsp, const cu_integer nCells,
    cu_real *__restrict__ hostDimMuKappa, cu2DArray<cu_real> &dDimMuKappa,
    const cudaStreams &streams, const bool copyDYiT) {
    if (copyDYiT) {
        dYiTP.copyFromHostAsync(hostYiTP, streams());
    }
    if (nTB.blockSize().x == 4) {
        getMuKappaDim<4><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(),
                           streams()>>>(dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else if (nTB.blockSize().x == 8) {
        getMuKappaDim<8><<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(),
                           streams()>>>(dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else if (nTB.blockSize().x == 16) {
        getMuKappaDim<16>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(),
               streams()>>>(dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else if (nTB.blockSize().x == 32) {
        getMuKappaDim<32>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(),
               streams()>>>(dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else if (nTB.blockSize().x == 64) {
        getMuKappaDim<64>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(),
               streams()>>>(dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else if (nTB.blockSize().x == 128) {
        getMuKappaDim<128>
            <<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(),
               streams()>>>(dYiTP, tran, dDimMuKappa, nsp, nCells);
    } else {
        fprintf(stderr,
                "OpenHurricane Error Of Calling CUDA function: Unsupported block "
                "size in x direction: %d at %s:%d of \"%s\" \n",
                nTB.blockSize().x, HUR_FILE, HUR_LINE, HUR_FUNCTION);
        exit(EXIT_FAILURE);
    }
    dDimMuKappa.copyToHostAsync(hostDimMuKappa, streams());
}

} // End of namespace CUDATransport
} // End of namespace OpenHurricane

#endif // CUDA_PARALLEL