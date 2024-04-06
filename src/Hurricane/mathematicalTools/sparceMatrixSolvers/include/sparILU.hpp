/*!
 * \file sparILU.hpp
 * \brief Header of sparse matrices Incomplete Lower Upper smoother for \f${\bf{A}}x = b\f$
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
         * \brief Incomplete Lower Upper solver for \f${\bf{A}}x = b\f$.
         * The Incomplete Lower Upper procedure can be illustrated as:
         * \f[\left( {{\bf{D}} + {\bf{L}}} \right){{\bf{D}}^{ - 1}}\left( {{\bf{D}} + {\bf{U}}} \right)\left( {{x^{\left( {k + 1} \right)}} - {x^{\left( k \right)}}} \right) = b - {\bf{A}}{x^{\left( k \right)}}\f]
         */
        template <class Type> class ILU : public sparSolver<Type> {
        public:
            using Base = sparSolver<Type>;
            using solution_type = typename Base::solution_type;

        private:
            /** \brief Inverse diagonal matrix. */
            mutable Array<Type> rD_;

            /** \brief D-DA */
            mutable Array<Type> DDA_;

        public:
            declareClassName(ILU);

            /**
             * \brief Construct from matrix.
             */
            inline ILU(const CRSMatrix<Type> &mat, const controller &cont)
                : Base(mat, cont), rD_(mat.nCells()), DDA_(mat.nCells()) {}

            /**
             * \brief Destructor.
             */
            virtual ~ILU() noexcept {}

        protected:
            void calcReciprocalD(Array<Type> &D, Array<Type> &rD, const CRSMatrix<Type> &A) const {
                const auto nCells = A.nCells();
                const auto &rowP = A.rowPtr();
                const auto &colId = A.col();
                const auto vv = A.data();

                const auto &cellNeberFace = A.addressing().cellNeberFace();

                const auto &lci = A.leftCellIndex();
                const auto &rci = A.rightCellIndex();

                for (integer i = 0; i < nCells; ++i) {
                    D[i] = A.Aii(i);
                    for (integer iNebr = 0; iNebr < cellNeberFace[i].size(); iNebr++) {
                        const auto fi = cellNeberFace[i][iNebr];
                        const auto II = colId[lci[fi]];
                        const auto JJ = colId[rci[fi]];

                        const auto &aJI = vv[lci[fi]];
                        const auto &aIJ = vv[rci[fi]];

                        const integer j = (i == II) ? JJ : II;
                        if (j < i) {
                            if (i == II) {
                                D[i] -= rD[j] * (aIJ * aJI);
                            } else // i == JJ
                            {
                                D[i] -= rD[j] * (aJI * aIJ);
                            }
                        }
                    }
                    rD[i] = inv(D[i]);
                    D[i] -= A.Aii(i);
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

                Array<solution_type> xTmp(nCells, Zero);

                solution_type bi;

                for (integer n = 0; n < nIter; ++n) {
                    Array<solution_type> bTmp(b);
                    A.interfaceTransfer(x, bTmp);

                    integer uEnd = rowP[0];
                    integer uStart;

                    // forward sweep
                    for (integer i = 0; i < nCells; ++i) {
                        // Start and end for this row of matrix A
                        uStart = uEnd;
                        uEnd = rowP[i + 1];
                        bTmp[i] += DDA_[i] * x[i];
                        bi = bTmp[i];

                        for (integer j0 = uStart; j0 < uEnd; ++j0) {
                            const auto &aij = vv[j0];
                            const auto j = colId[j0];
                            if (j < i) // Lower tridiagonal parts for this row of matrix A: j < i
                            {
                                bi -= aij * xTmp[j];
                            } else if (j >
                                       i) // Upper tridiagonal parts for this row of matrix A: j > i
                            {
                                bi -= aij * x[j];
                            }
                        }

                        xTmp[i] = rD_[i] * bi;
                    }

                    uStart = rowP[nCells];
                    // backward sweep
                    for (integer i = nCells - 1; i >= 0; --i) {
                        // Start and end for this row of matrix A
                        uEnd = uStart;
                        uStart = rowP[i];

                        bi = bTmp[i];

                        for (integer j0 = uStart; j0 < uEnd; ++j0) {
                            const auto &aij = vv[j0];
                            const auto j = colId[j0];
                            if (j < i) // Lower tridiagonal parts for this row of matrix A: j < i
                            {
                                bi -= aij * xTmp[j];
                            } else if (j >
                                       i) // Upper tridiagonal parts for this row of matrix A: j > i
                            {
                                bi -= aij * x[j];
                            }
                        }
                        x[i] = rD_[i] * bi;
                    }
                }
            }

        public:
            virtual inline void setupSolver() const { calcReciprocalD(DDA_, rD_, Base::A_); }

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
        using realILUSolver = ILU<real>;
        using realBlockILUSolver = ILU<realSquareMatrix>;
    } // namespace sparseSolver
} // namespace OpenHurricane
