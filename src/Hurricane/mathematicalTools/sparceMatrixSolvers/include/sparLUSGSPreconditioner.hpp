/*!
 * \file sparLUSGSPreconditioner.hpp
 * \brief Header of LUSGS Preconditioner.
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
#include "sparPreconditioner.hpp"

namespace OpenHurricane {
    namespace sparseSolver {

        template <class Type> class LUSGSPreconditioner : public preconditioner<Type> {
        public:
            using Base = preconditioner<Type>;

        private:
            mutable Array<Type> invA_;

            void calcReciprocal() const {
                const auto &A = Base::A_;
                bool hasInversedDiag = A.hasInversedDiag();
                if (!hasInversedDiag) {
                    const auto nCells = A.nCells();
                    invA_.resize(nCells);
                    for (integer i = 0; i < nCells; ++i) {
                        const auto &aii = A.Aii(i);
                        invA_[i] = inv(aii);
                    }
                }
            }

        public:
            using value_type = typename Base::value_type;
            using solution_type = typename Base::solution_type;

            declareClassName(LUSGSPreconditioner);

            inline LUSGSPreconditioner(const CRSMatrix<Type> &mat, const controller &cont)
                : Base(mat, cont), invA_() {
                calcReciprocal();
            }

            inline virtual ~LUSGSPreconditioner() noexcept {}

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             */
            virtual void precondition(Array<solution_type> &x,
                                      const Array<solution_type> &x0) const {
                const auto &A = Base::A_;
                const auto nCells = A.nCells();
                const auto &rowP = A.rowPtr();
                const auto &colId = A.col();
                const auto vv = A.data();

                solution_type bi;

                bool hasInversedDiag = A.hasInversedDiag();

                integer uEnd = rowP[0];
                integer uStart;

                for (integer i = 0; i < nCells; ++i) {
                    uStart = uEnd;
                    uEnd = rowP[i + 1];
                    bi = x0[i];

                    for (integer j = uStart; j < uEnd; ++j) {
                        if (colId[j] < i) {
                            bi -= vv[j] * x[colId[j]];
                        }
                    }

                    if (hasInversedDiag) {
                        x[i] = A.invAii(i) * bi;
                    } else {
                        x[i] = invA_[i] * bi;
                    }
                }

                uStart = rowP[nCells];
                for (integer i = nCells - 1; i >= 0; --i) {
                    uEnd = uStart;
                    uStart = rowP[i];

                    bi = Zero;
                    for (integer j = uStart; j < uEnd; ++j) {
                        if (colId[j] > i) {
                            bi += vv[j] * x[colId[j]];
                        }
                    }

                    if (hasInversedDiag) {
                        x[i] -= A.invAii(i) * bi;
                    } else {
                        x[i] -= invA_[i] * bi;
                    }
                }
            }
        };

        using realLUSGSPreconditioner = LUSGSPreconditioner<real>;
        using realBlockLUSGSPreconditioner = LUSGSPreconditioner<realSquareMatrix>;

    } // namespace sparseSolver
} // namespace OpenHurricane