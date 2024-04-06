/*!
 * \file cuSparJacobi.cpp
 * \brief Main subroutines for sparse matrices Jacobi solver for \f${\bf{A}}x = b\f$ on CUDA platform.
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
#include "cuSparJacobi.hpp"
#ifdef CUDA_PARALLEL

#include "cuSparJacobi.h"
#include "nThreadsAndBlocks.hpp"

namespace OpenHurricane {
    namespace sparseSolver {
        template <class CSRSMatType>
        createClassNameStr(cuSparJacobiBlockSolver<CSRSMatType>, "Jacobi");

        createClassNameStrTmpl(cuRealMBlockSparJacobiSolver, "Jacobi");

        registerObjFty(cuRealMBlockSparSolver, cuRealMBlockSparJacobiSolver, controller);
        registerObjFty(cuRealMBlockSparSolver, cuRealMBlockSparJacobiSolver, Smoother);

        template <>
        void cuRealMBlockSparJacobiSolver::smooth(cuMBlockArray<cu_real, cu_integer> &x,
                                                 const cuMCRSBlockMatrix<cu_real> &A,
                                                 const cuMBlockArray<cu_real, cu_integer> &b,
                                                 const integer nIter, const int devId,
                                                 cudaStream_t s) const {
            OpenHurricane::nThreadsAndBlocks nTB(x.blockDimy(), A.nRows());
            auto &nB = nTB.blockSize();
            nTB.setSharedMemSize1<cu_real>(nB.y * x.blockDimy());
            JacobiaMBlockSmooth(x, A, b, nIter, nTB, devId, s);
        }

        createClassNameStrTmpl(cuRealBlockSparJacobiSolver, "Jacobi");

        registerObjFty(cuRealBlockSparSolver, cuRealBlockSparJacobiSolver, controller);
        registerObjFty(cuRealBlockSparSolver, cuRealBlockSparJacobiSolver, Smoother);

        template <>
        void cuRealBlockSparJacobiSolver::smooth(cuBlockArray1<cu_real, cu_integer> &x,
                                                  const cuCRSBlockMatrix<cu_real> &A,
                                                  const cuBlockArray1<cu_real, cu_integer> &b,
                                                  const integer nIter, const int devId,
                                                  cudaStream_t s) const {
            OpenHurricane::nThreadsAndBlocks nTB(x.blockDimy(), A.nRows());
            nTB.setBlockAndGridSize(x.blockDimy(), A.nRows());
            auto &nB = nTB.blockSize();
            nTB.setSharedMemSize1<cu_real>(nB.y * x.blockDimy());
            JacobiaBlockSmooth(x, A, b, nIter, nTB, devId, s);
        }
    } // namespace sparseSolver
} // namespace OpenHurricane
#endif //CUDA_PARALLEL
