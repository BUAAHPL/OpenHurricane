/*!
 * \file cuSparJacobi.hpp
 * \brief Header of sparse matrices Jacobi solver for \f${\bf{A}}x = b\f$ on CUDA platform
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
#include "cuSparSolver.hpp"

#ifdef CUDA_PARALLEL

namespace OpenHurricane {

    namespace sparseSolver {
        /**
         * \brief Abstract base-class for sparse block matrix solver for  \f${\bf{A}}x = b\f$ on CUDA platform.
         */
        template <class CSRSMatType>
        class cuSparJacobiBlockSolver : public cuSparBlockSolver<CSRSMatType> {
        public:
            using Base = cuSparBlockSolver<CSRSMatType>;
            using value_type = typename Base::value_type;
            using solution_type = typename Base::solution_type;

        public:
            declareClassNames;

            /**
             * \brief Construct from matrix.
             */
            inline cuSparJacobiBlockSolver(const CSRSMatType &mat,
                                           const controller &cont)
                : Base(mat, cont) {}

            /**
             * \brief Destructor.
             */
            virtual ~cuSparJacobiBlockSolver() noexcept {}

        protected:
            void calcDiagA(cudaStream_t s = 0) const {
                //LFatal("This function has not been implemented");
            }
            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             */
            void smooth(solution_type &x, const CSRSMatType &A, const solution_type &b,
                        const cu_integer nIter,
                        const int devId, cudaStream_t s = 0) const {
                LFatal("This function has not been implemented");
            }

        public:
            virtual void setupSolver(cudaStream_t s = 0) const {}

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             */
            virtual void solve(solution_type &x, const solution_type &b, const cu_integer nIter,
                               const int devId, cudaStream_t s = 0) const {
                setupSolver(s);
                smooth(x, Base::A_, b, nIter, devId, s);
            }

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             * Must be called after setupSolver().
             */
            virtual void solveOnly(solution_type &x, const solution_type &b, const cu_integer nIter,
                                   const int devId, cudaStream_t s = 0) const {
                smooth(x, Base::A_, b, nIter, devId, s);
            }
        };

        using cuRealMBlockSparJacobiSolver = cuSparJacobiBlockSolver<cuMCRSBlockMatrix<cu_real>>;
        using cuRealBlockSparJacobiSolver = cuSparJacobiBlockSolver<cuCRSBlockMatrix<cu_real>>;

        template <>
        void cuRealMBlockSparJacobiSolver::smooth(cuMBlockArray<cu_real, cu_integer> &x,
                                                 const cuMCRSBlockMatrix<cu_real> &A,
                                                 const cuMBlockArray<cu_real, cu_integer> &b,
                                                 const integer nIter, const int devId,
                                                  cudaStream_t s) const;
        template <>
        void cuRealBlockSparJacobiSolver::smooth(cuBlockArray1<cu_real, cu_integer> &x,
                                                  const cuCRSBlockMatrix<cu_real> &A,
                                                  const cuBlockArray1<cu_real, cu_integer> &b,
                                                  const integer nIter, const int devId,
                                                  cudaStream_t s) const;

    } // namespace sparseSolver
} // namespace OpenHurricane
#endif // CUDA_PARALLEL