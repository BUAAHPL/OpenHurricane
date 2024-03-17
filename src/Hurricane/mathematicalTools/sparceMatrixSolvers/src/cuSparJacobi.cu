/*!
 * \file cuSparJacobi.cu
 * \brief Main subroutines for sparse matrices Jacobi solver for \f${\bf{A}}x = b\f$ on CUDA platform.
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
#ifdef CUDA_PARALLEL
#include "cuSparJacobi.h"

namespace OpenHurricane {
    namespace sparseSolver {

        cu_global void JacobiMBlockSolver(const cuMBlockArray<cu_real, cu_integer> x,
                                          cuMBlockArray<cu_real, cu_integer> newX,
                                          const cuMCRSBlockMatrix<cu_real> A,
                                          const cuMBlockArray<cu_real, cu_integer> b) {
            // The index of this cell
            const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;
            const auto nBlocky = x.blockDimy();

            const auto nCells = A.nRows();
            const auto &rowP = A.rowPtr();
            const auto &colId = A.col();
            const auto &vv = A.value();
            const auto &diagId = A.diagIndex();

            extern __shared__ volatile char schar[];
            volatile cu_real *bis;
            bis = (cu_real *)&schar[0];

            volatile cu_real *bi = (cu_real *)(bis + threadIdx.y * nBlocky);

            auto tid = threadIdx.x;
            if (cellId < nCells) {
                const auto b0 = b(cellId);
                while (tid < nBlocky) {
                    bi[tid] = b0(0, tid);
                    tid += blockDim.x;
                }
            }

            __syncthreads();

            if (cellId < nCells) {
                cu_integer uEnd = rowP[cellId + 1];
                cu_integer uStart = rowP[cellId];
                for (cu_integer j0 = uStart; j0 < uEnd; ++j0) {
                    const auto j = colId[j0];
                    if (j != cellId) {
                        const auto &aij = vv(j0);
                        const auto x0 = x(j);
                        tid = threadIdx.x;
                        while (tid < nBlocky) {
                            for (cu_ushort kk = 0; kk < nBlocky; ++kk) {
                                bi[tid] -= aij(tid, kk) * x0(0, kk);
                            }
                            tid += blockDim.x;
                        }
                    }
                }
            }
            __syncthreads();

            if (cellId < nCells) {
                auto xi = newX(cellId);
                const auto &aii = vv(diagId[cellId]);
                tid = threadIdx.x;
                while (tid < nBlocky) {
                    xi(0, tid) = bi[tid] / aii(tid, tid);
                    tid += blockDim.x;
                }
            }
        }

        void JacobiaMBlockSmooth(cuMBlockArray<cu_real, cu_integer> &x,
                                 const cuMCRSBlockMatrix<cu_real> &A,
                                 const cuMBlockArray<cu_real, cu_integer> &b,
                                 const cu_integer nIter, nThreadsAndBlocks &nTB, const int devId,
                                 cudaStream_t s) {
            cuMBlockArray<cu_real, cu_integer> newX(x.size(), 1, x.blockDimy());

            for (cu_integer n = 0; n < nIter; ++n) {
                JacobiMBlockSolver<<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), s>>>(
                    x, newX, A, b);
                newX.memPrefetchAsync(cudaCpuDeviceId, s);
                checkCUDAError(cudaStreamSynchronize(s));
                for (cu_integer i = 0; i < x.size(); ++i) {
                    for (cu_ushort j = 0; j < x.blockDimy(); ++j) {
                        *(x(i).data() + j) = *(newX(i).data() + j);
                    }
                }
                if (n != nIter - 1) {
                    x.memPrefetchAsync(devId, s);
                }
            }
            destroyCuArray(newX);
        }

        cu_global void JacobiBlockSolver(const cuBlockArray1<cu_real, cu_integer> x,
                                         cuBlockArray1<cu_real, cu_integer> newX,
                                         const cuCRSBlockMatrix<cu_real> A,
                                         const cuBlockArray1<cu_real, cu_integer> b) {
            // The index of this cell
            const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;
            const auto nBlocky = x.blockDimy();

            const auto nCells = A.nRows();
            const auto &rowP = A.rowPtr();
            const auto &colId = A.col();
            const auto &vv = A.value();
            const auto &diagId = A.diagIndex();

            extern __shared__ volatile char schar[];
            volatile cu_real *bis;
            bis = (cu_real *)&schar[0];

            volatile cu_real *bi = (cu_real *)(bis + threadIdx.y * nBlocky);

            auto tid = threadIdx.x;
            if (cellId < nCells) {
                const auto b0 = b(cellId);
                while (tid < nBlocky) {
                    bi[tid] = b0(0, tid);
                    tid += blockDim.x;
                }
            }

            __syncthreads();

            if (cellId < nCells) {
                cu_integer uEnd = rowP(cellId + 1);
                cu_integer uStart = rowP(cellId);

                for (cu_integer j0 = uStart; j0 < uEnd; ++j0) {
                    const auto j = colId(j0);
                    if (j != cellId) {
                        const auto &aij = vv(j0);
                        const auto x0 = x(j);
                        tid = threadIdx.x;
                        while (tid < nBlocky) {
                            for (cu_ushort kk = 0; kk < nBlocky; ++kk) {
                                bi[tid] -= aij(tid, kk) * x0(0, kk);
                            }
                            tid += blockDim.x;
                        }
                    }
                }
            }
            __syncthreads();

            if (cellId < nCells) {
                auto xi = newX(cellId);
                const auto &aii = vv(diagId(cellId));
                tid = threadIdx.x;
                while (tid < nBlocky) {
                    xi(0, tid) = bi[tid] / aii(tid, tid);
                    tid += blockDim.x;
                }
            }  
        }

        cu_global void getNewSolution(cuBlockArray1<cu_real, cu_integer> x,
                                      const cuBlockArray1<cu_real, cu_integer> newX) {
            // The index of this cell
            const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;
            const auto nBlocky = x.blockDimy();
            if (cellId < x.size()) {
                auto xi = newX(cellId);
                auto xn = x(cellId);
                auto tid = threadIdx.x;
                while (tid < nBlocky) {
                    xn(0, tid) = xi(0, tid);
                    tid += blockDim.x;
                }
            }
        }

        cu_global void setZeroSolution(cuBlockArray1<cu_real, cu_integer> x) {
            // The index of this cell
            const unsigned int cellId = threadIdx.y + blockDim.y * blockIdx.x;
            const auto nBlocky = x.blockDimy();
            if (cellId < x.size()) {
                auto xn = x(cellId);
                auto tid = threadIdx.x;
                while (tid < nBlocky) {
                    xn(0, tid) = 0;
                    tid += blockDim.x;
                }
            }
        }

        void JacobiaBlockSmooth(cuBlockArray1<cu_real, cu_integer> &x,
                                const cuCRSBlockMatrix<cu_real> &A,
                                const cuBlockArray1<cu_real, cu_integer> &b, const cu_integer nIter,
                                nThreadsAndBlocks &nTB, const int devId, cudaStream_t s) {
            setZeroSolution<<<nTB.gridSize(), nTB.blockSize(), 0, s>>>(x);

            cuBlockArray1<cu_real, cu_integer> newX(x.size(), 1, x.blockDimy());

            for (cu_integer n = 0; n < nIter; ++n) {
                JacobiBlockSolver<<<nTB.gridSize(), nTB.blockSize(), nTB.sharedMemSize(), s>>>(
                    x, newX, A, b);
                getNewSolution<<<nTB.gridSize(), nTB.blockSize(), 0, s>>>(x, newX);
            }

            x.copyToHostAsync(s);

            checkCUDAError(cudaStreamSynchronize(s));
            destroyCuArray(newX);
        }
    } // namespace sparseSolver
} // namespace OpenHurricane
#endif //CUDA_PARALLEL
