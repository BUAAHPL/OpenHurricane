/*!
 * \file sparGaussSeidel.hpp
 * \brief Header of sparse matrices Gauss-Seidel Smoother for \f${\bf{A}}x = b\f$
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
#include "sparSolver.hpp"

namespace OpenHurricane {
    namespace sparseSolver {
        /**
         * \brief Gauss-Seidel solver for \f${\bf{A}}x = b\f$.
         * The Gauss-Seidel procedure can be illustrated using the scalar system, the sweep can be written as:
         * \f[x_i^{\left( {k + 1} \right)} = {{\left( {{b_i} - \sum\limits_{j < i} {{a_{ij}}x_j^{\left( {k + 1} \right)}}  - \sum\limits_{j > i} {{a_{ij}}x_j^{\left( k \right)}} } \right)} \mathord{\left/
                {\vphantom {{\left( {{b_i} - \sum\limits_{j < i} {{a_{ij}}x_j^{\left( {k + 1} \right)}}  - \sum\limits_{j > i} {{a_{ij}}x_j^{\left( k \right)}} } \right)} {{a_{ii}}}}} \right.
                \kern-\nulldelimiterspace} {{a_{ii}}}},{\rm{  }}\left( {i = 1,2, \ldots ,n;k = 0,1, \ldots } \right)\f]
         */
        template <class Type> class GaussSeidel : public sparSolver<Type> {
        public:
            using Base = sparSolver<Type>;
            using solution_type = typename Base::solution_type;

        public:
            declareClassName(GaussSeidel);

            /**
             * \brief Construct from matrix.
             */
            inline GaussSeidel(const CRSMatrix<Type> &mat, const controller &cont)
                : Base(mat, cont), invA_() {}

            /**
             * \brief Destructor.
             */
            virtual ~GaussSeidel() noexcept {}

        protected:
            mutable List<Type> invA_;

            void calcDiagA() const {
                bool hasInversedDiag = Base::A_.hasInversedDiag();
                if (hasInversedDiag) {
                    return;
                }
                const auto nCells = Base::A_.nCells();
                if (!hasInversedDiag) {
                    invA_.resize(nCells);
                }
                for (integer i = 0; i < nCells; ++i) {
                    if (!hasInversedDiag) {
                        const Type &aii = Base::A_.Aii(i);
                        invA_[i] = inv(aii);
                    }
                }
            }

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             */
            void smooth(Array<solution_type> &x, const CRSMatrix<Type> &A,
                        const Array<solution_type> &b, const integer nIter) const {
                const auto nCells = A.nCells();
                const auto &rowP = A.rowPtr();
                const auto &colId = A.col();
                const auto vv = A.data();

                solution_type bi;

                bool hasInversedDiag = A.hasInversedDiag();

                for (integer n = 0; n < nIter; ++n) {
                    Array<solution_type> bTmp(b);
                    A.interfaceTransfer(x, bTmp);
                    integer uEnd = rowP[0];
                    integer uStart;

                    for (integer i = 0; i < nCells; ++i) {
                        // Start and end for this row of matrix A
                        uStart = uEnd;
                        uEnd = rowP[i + 1];

                        bi = bTmp[i];
                        for (integer j0 = uStart; j0 < uEnd; ++j0) {
                            const auto &aij = vv[j0];
                            const auto j = colId[j0];
                            if (j < i) // Lower tridiagonal parts for this row of matrix A: j < i
                            {
                                // x_j^(k+1)
                                bi -= aij * x[j];
                            } else if (j >
                                       i) // Upper tridiagonal parts for this row of matrix A: j > i
                            {
                                // x_j^(k)
                                bi -= aij * x[j];
                            }
                        }

                        // x_i^(k+1) = bi/Aii
                        if (hasInversedDiag) {
                            x[i] = A.invAii(i) * bi;
                        } else {
                            x[i] = invA_[i] * bi;
                        }
                    }
                }
            }

        public:
            virtual inline void setupSolver() const { calcDiagA(); }

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             */
            virtual void solve(Array<solution_type> &x, const Array<solution_type> &b,
                               const integer nIter) const {
                setupSolver();
                smooth(x, Base::A_, b, nIter);
            }

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             * Must be called after setupSolver().
             */
            virtual void solveOnly(Array<solution_type> &x, const Array<solution_type> &b,
                                   const integer nIter) const {
                smooth(x, Base::A_, b, nIter);
            }
        };

        using realGaussSeidelSolver = GaussSeidel<real>;
        using realBlockGaussSeidelSolver = GaussSeidel<realSquareMatrix>;
    } // namespace sparseSolver
} // namespace OpenHurricane
