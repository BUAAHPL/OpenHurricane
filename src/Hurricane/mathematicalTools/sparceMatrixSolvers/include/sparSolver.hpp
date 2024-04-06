/*!
 * \file sparSolver.hpp
 * \brief Header of sparse matrices solver for \f${\bf{A}}x = b\f$
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
#include "CRSMatrix.hpp"

#include "controller.hpp"
#include "objectFactory.hpp"

namespace OpenHurricane {
    namespace sparseSolver {
        /**
         * \brief Abstract base-class for sparse matrix smoother for  \f${\bf{A}}x = b\f$.
         */
        template <class Type> class sparSolver {
        public:
            using value_type = typename CRSMatrix<Type>::value_type;
            using solution_type = typename CRSMatrix<Type>::solution_type;

        protected:
            const CRSMatrix<Type> &A_;

            real RTol_;
            real ATol_;

        public:
            declareClassNames;
            declareObjFty(sparSolver, controller,
                          (const CRSMatrix<Type> &mat, const controller &cont), (mat, cont));

            declareObjFty(sparSolver, Smoother,
                          (const CRSMatrix<Type> &mat, const controller &cont), (mat, cont));

            /**
             * \brief Construct from matrix.
             */
            inline sparSolver(const CRSMatrix<Type> &mat, const controller &cont)
                : A_(mat), RTol_(cont.findOrDefault<real>("rtol", 1e-4)),
                  ATol_(cont.findOrDefault<real>("atol", 1e-15)) {}

            hur_nodiscard static uniquePtr<sparSolver<Type>> creator(const CRSMatrix<Type> &mat,
                                                                     const controller &cont) {
                string solverType = cont.findWord(sparSolver<Type>::className_);
                defineInObjCreatorTmpl(sparSolver<Type>, solverType, controller, (mat, cont));
            }
            
            hur_nodiscard static uniquePtr<sparSolver<Type>> createSmoother(const CRSMatrix<Type> &mat,
                                                                     const controller &cont) {
                string solverType = cont.findWord(sparSolver<Type>::className_);
                defineInObjCreatorTmpl(sparSolver<Type>, solverType, Smoother, (mat, cont));
            }


            /**
             * \brief Destructor.
             */
            virtual ~sparSolver() noexcept {}

            virtual void setupSolver() const {}

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             */
            virtual void solve(Array<solution_type> &x, const Array<solution_type> &b,
                               const integer nIter) const = 0;
            
            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             * Must be called after setupSolver().
             */
            virtual void solveOnly(Array<solution_type> &x, const Array<solution_type> &b,
                               const integer nIter) const = 0;
        };

        using realSparSolver = sparSolver<real>;
        using realBlockSparSolver = sparSolver<realSquareMatrix>;
    } // namespace sparseSolver
} // namespace OpenHurricane
