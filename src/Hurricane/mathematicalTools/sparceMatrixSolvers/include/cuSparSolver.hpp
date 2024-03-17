/*!
 * \file cuSparSolver.hpp
 * \brief Header of sparse matrices solver for \f${\bf{A}}x = b\f$ on CUDA platform
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
#include "cuArrayBase.hpp"
#include "cuCRSMatrix.hpp"

#include "controller.hpp"
#include "objectFactory.hpp"

#ifdef CUDA_PARALLEL
#include "cudaStreams.hpp"

namespace OpenHurricane {

    namespace sparseSolver {
        /**
         * \brief Abstract base-class for sparse block matrix solver for  \f${\bf{A}}x = b\f$ on CUDA platform.
         */
        template <class CSRSMatType> class cuSparBlockSolver {
        public:
            using value_type = typename CSRSMatType::value_type;
            using solution_type = typename CSRSMatType::solution_type;

        protected:
            const CSRSMatType &A_;

            real RTol_;
            real ATol_;

        public:
            declareClassNames;
            declareObjFty(cuSparBlockSolver, controller,
                          (const CSRSMatType &mat, const controller &cont), (mat, cont));

            declareObjFty(cuSparBlockSolver, Smoother,
                          (const CSRSMatType &mat, const controller &cont), (mat, cont));

            /**
             * \brief Construct from matrix.
             */
            inline cuSparBlockSolver(const CSRSMatType &mat, const controller &cont)
                : A_(mat), RTol_(cont.findOrDefault<real>("rtol", 1e-4)),
                  ATol_(cont.findOrDefault<real>("atol", 1e-15)) {}

            hur_nodiscard static uniquePtr<cuSparBlockSolver<CSRSMatType>>
            creator(const CSRSMatType &mat, const controller &cont) {
                string solverType = cont.findWord(cuSparBlockSolver<CSRSMatType>::className_);
                defineInObjCreatorTmpl(cuSparBlockSolver<CSRSMatType>, solverType, controller,
                                       (mat, cont));
            }

            hur_nodiscard static uniquePtr<cuSparBlockSolver<CSRSMatType>>
            createSmoother(const CSRSMatType &mat, const controller &cont) {
                string solverType = cont.findWord(cuSparBlockSolver<CSRSMatType>::className_);
                defineInObjCreatorTmpl(cuSparBlockSolver<CSRSMatType>, solverType, Smoother,
                                       (mat, cont));
            }

            /**
             * \brief Destructor.
             */
            virtual ~cuSparBlockSolver() noexcept {}

            virtual void setupSolver(cudaStream_t s = 0) const {}

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             */
            virtual void solve(solution_type &x, const solution_type &b, const cu_integer nIter,
                               const int devId, cudaStream_t s = 0) const = 0;

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             * Must be called after setupSolver().
             */
            virtual void solveOnly(solution_type &x, const solution_type &b, const cu_integer nIter,
                                   const int devId, cudaStream_t s = 0) const = 0;
        };

        using cuRealMBlockSparSolver = cuSparBlockSolver<cuMCRSBlockMatrix<cu_real>>;
        using cuRealBlockSparSolver = cuSparBlockSolver<cuCRSBlockMatrix<cu_real>>;
    } // namespace sparseSolver
} // namespace OpenHurricane
#endif // CUDA_PARALLEL