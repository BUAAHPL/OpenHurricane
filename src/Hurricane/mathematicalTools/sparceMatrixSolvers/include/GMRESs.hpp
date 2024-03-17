/*!
 * \file GMRESs.hpp
 * \brief Header of spar GMRES solver
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
#include "sparPreconditioner.hpp"
#include "sparSolver.hpp"

namespace OpenHurricane {

    namespace sparseSolver {
        template <class Type> class GMRES : public sparSolver<Type> {
        public:
            using solution_type = typename sparSolver<Type>::solution_type;
            using Base = sparSolver<Type>;

        private:
            const controller &gmCont_;
            uniquePtr<sparseSolver::sparSolver<Type>> smootherPtr_;
            integer nFirstSolve_;
            integer krylovDimension_;

            real epsilon_;

            integer nRestarts_;

            hur_nodiscard solution_type L2Norm(const Array<solution_type> &x) const {
                LFatal("Attempt to call null function");
            }

            hur_nodiscard solution_type dotPro(const Array<solution_type> &x,
                                               const Array<solution_type> &y) const {
                LFatal("Attempt to call null function");
            }

            void update(const Array<Array<solution_type>> &H, Array<solution_type> &x,
                        Array<solution_type> &zets, const Array<Array<solution_type>> &v,
                        integer iDim) const {
                LFatal("Attempt to call null function");
            }

        public:
            declareClassName(GMRES);

            inline GMRES(const CRSMatrix<Type> &mat, const controller &cont)
                : Base(mat, cont), gmCont_(cont.subController("GMRES")), smootherPtr_(nullptr),
                  nFirstSolve_(2), krylovDimension_(40), epsilon_(1e-7), nRestarts_(3) {
                if (gmCont_.found("sparSmoother")) {
                    smootherPtr_ = sparseSolver::sparSolver<Type>::createSmoother(mat, gmCont_);
                    nFirstSolve_ = gmCont_.findOrDefault<integer>("nFirstSolve", nFirstSolve_);
                }
                krylovDimension_ =
                    gmCont_.findOrDefault<integer>("krylovDimension", krylovDimension_);
                epsilon_ = gmCont_.findOrDefault<real>("epsilon", epsilon_);
                nRestarts_ = gmCont_.findOrDefault<integer>("nRestarts", nRestarts_);
            }

            inline virtual ~GMRES() noexcept { smootherPtr_.clear(); }
                        
            virtual void setupSolver() const {}

            virtual void solve(Array<solution_type> &x, const Array<solution_type> &b,
                               const integer nIter) const {
                LFatal("Attempt to call null function");
            }

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             * Must be called after setupSolver().
             */
            virtual void solveOnly(Array<solution_type> &x, const Array<solution_type> &b,
                                   const integer nIter) const {
                LFatal("Attempt to call null function");
            }
        };

        using sparRealGMRESSolver = GMRES<real>;
        using sparRealBlockGMRESSolver = GMRES<realSquareMatrix>;

        template <>
        hur_nodiscard GMRES<real>::solution_type
        GMRES<real>::L2Norm(const Array<solution_type> &x) const;

        template <>
        hur_nodiscard GMRES<real>::solution_type
        GMRES<real>::dotPro(const Array<solution_type> &x, const Array<solution_type> &y) const;

        template <>
        void GMRES<real>::update(const Array<Array<solution_type>> &H, Array<solution_type> &x,
                                 Array<solution_type> &zets, const Array<Array<solution_type>> &v,
                                 integer iDim) const;

        template <>
        void GMRES<real>::solve(Array<solution_type> &x, const Array<solution_type> &b,
                                const integer nIter) const;

        template <>
        hur_nodiscard GMRES<realSquareMatrix>::solution_type
        GMRES<realSquareMatrix>::L2Norm(const Array<solution_type> &x) const;

        template <>
        hur_nodiscard GMRES<realSquareMatrix>::solution_type
        GMRES<realSquareMatrix>::dotPro(const Array<solution_type> &x,
                                        const Array<solution_type> &y) const;

        template <>
        void GMRES<realSquareMatrix>::update(const Array<Array<solution_type>> &H,
                                             Array<solution_type> &x, Array<solution_type> &zets,
                                             const Array<Array<solution_type>> &v,
                                             integer iDim) const;

        template <>
        void GMRES<realSquareMatrix>::solve(Array<solution_type> &x, const Array<solution_type> &b,
                                            const integer nIter) const;
    } // namespace sparseSolver
} // namespace OpenHurricane