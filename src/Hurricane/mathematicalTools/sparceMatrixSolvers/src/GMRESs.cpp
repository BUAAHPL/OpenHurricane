/*!
 * \file GMRESs.cpp
 * \brief Main subroutines for GMRES solver
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
#include "GMRESs.hpp"

namespace OpenHurricane {
    namespace sparseSolver {
        createClassNameStrTmpl(sparRealGMRESSolver, "GMRES");

        createClassNameStrTmpl(sparRealBlockGMRESSolver, "GMRES");

        registerObjFty(realSparSolver, sparRealGMRESSolver, controller);
        registerObjFty(realBlockSparSolver, sparRealBlockGMRESSolver, controller);

        template <>
        hur_nodiscard GMRES<real>::solution_type
        GMRES<real>::L2Norm(const Array<solution_type> &x) const {
            solution_type tmp = 0;

            for (const auto &e : x) {
                tmp += sqr(e);
            }
            return sqrt(tmp);
        }

        template <>
        hur_nodiscard GMRES<real>::solution_type
        GMRES<real>::dotPro(const Array<solution_type> &x, const Array<solution_type> &y) const {
            solution_type tmp = 0;
            for (integer i = 0; i < A_.nCells(); ++i) {
                tmp += x[i] * y[i];
            }
            return tmp;
        }

        template <>
        void GMRES<real>::update(const Array<Array<solution_type>> &H, Array<solution_type> &x,
                                 Array<solution_type> &zets, const Array<Array<solution_type>> &v,
                                 integer iDim) const {
            for (integer i = iDim; i >= 0; --i) {
                zets[i] /= H[i][i];
                for (integer j = i - 1; j >= 0; --j) {
                    zets[j] -= zets[i] * H[j][i];
                }
            }
            for (integer i = 0; i <= iDim; ++i) {
                for (integer n = 0; n < A_.nCells(); ++n) {
                    x[n] += v[i][n] * zets[i];
                }
            }
        }

        template <>
        void GMRES<real>::solve(Array<solution_type> &x, const Array<solution_type> &b,
                                const integer nIter) const {
            if (smootherPtr_) {
                smootherPtr_->solve(x, b, nFirstSolve_);
            }
            auto precond = sparseSolver::preconditioner<real>::creator(A_, gmCont_);

            Array<solution_type> residual(A_.nCells());
            Array<solution_type> Ax(A_.nCells());
            A_.multiple(x, Ax);
            precond->precondition(residual, b - Ax);

            auto beta = L2Norm(residual);

            if (beta < epsilon_) {
                return;
            }

            // Allocate the space for the givens rotations, and the upper
            // Hessenburg matrix.
            Array<Array<solution_type>> H(krylovDimension_ + 1);
            for (auto &e : H) {
                e.resize(krylovDimension_);
            }

            // The Givens rotations include the sine and cosine term. The
            // cosine term is in column zero, and the sine term is in column
            // one.
            Array<vector2D> givens(krylovDimension_ + 1);

            Array<Array<solution_type>> v(krylovDimension_ + 1);
            for (auto &e : v) {
                e.resize(A_.nCells());
            }

            Array<solution_type> zeta(krylovDimension_ + 1);

            real relres = 0;
            for (integer k = 0; k < nRestarts_; ++k) {
                v[0] = residual * inv(beta);
                zeta.setZero();
                zeta[0] = beta;
                integer m = 0;
                for (integer j = 0; j < krylovDimension_; ++j) {
                    A_.multiple(v[j], Ax);
                    precond->precondition(v[j + 1], Ax);

                    m = j + 1;
                    for (integer i = 0; i <= j; ++i) {
                        H[i][j] = dotPro(v[i], v[j + 1]);
                        v[j + 1] -= (H[i][j] * v[i]);
                    }

                    H[j + 1][j] = L2Norm(v[j + 1]);
                    if (H[j + 1][j] == 0) {
                        break;
                    }
                    v[j + 1] *= inv(H[j + 1][j]);

                    for (integer i = 0; i < j; ++i) {
                        auto tmp = givens[i][0] * H[i][j] + givens[i][1] * H[i + 1][j];
                        H[i + 1][j] = -givens[i][1] * H[i][j] + givens[i][0] * H[i + 1][j];
                        H[i][j] = tmp;
                    }

                    if (H[j + 1][j] == 0) {
                        givens[j][0] = 1;
                        givens[j][1] = 0;
                    } else if (mag(H[j][j]) > mag(H[j + 1][j])) {
                        auto tau = H[j + 1][j] / H[j][j];
                        givens[j][0] = inv(sqrt(real(1) + sqr(tau)));
                        givens[j][1] = givens[j][0] * tau;
                    } else {
                        auto tau = H[j][j] / H[j + 1][j];
                        givens[j][1] = inv(sqrt(real(1) + sqr(tau)));
                        givens[j][0] = givens[j][1] * tau;
                    }

                    auto hjj = H[j][j];
                    auto hjp1j = H[j + 1][j];

                    H[j][j] = givens[j][0] * hjj + givens[j][1] * hjp1j;
                    H[j + 1][j] = -givens[j][1] * hjj + givens[j][0] * hjp1j;

                    auto zj = zeta[j];
                    auto zjp1 = zeta[j + 1];
                    zeta[j] = givens[j][0] * zj + givens[j][1] * zjp1;
                    zeta[j + 1] = -givens[j][1] * zj + givens[j][0] * zjp1;

                    relres = mag(zeta[j + 1]);
                    if (relres < epsilon_ * beta) {
                        break;
                    }
                }

                update(H, x, zeta, v, m - 1);

                A_.multiple(x, Ax);

                precond->precondition(residual, b - Ax);

                beta = L2Norm(residual);
            }
        }

        template <>
        hur_nodiscard GMRES<realSquareMatrix>::solution_type
        GMRES<realSquareMatrix>::L2Norm(const Array<solution_type> &x) const {
            solution_type tmp(x.first().size(), Zero);

            for (const auto &e : x) {
                tmp += sqr(e);
            }
            return sqrt(tmp);
        }

        template <>
        hur_nodiscard GMRES<realSquareMatrix>::solution_type
        GMRES<realSquareMatrix>::dotPro(const Array<solution_type> &x,
                                        const Array<solution_type> &y) const {
            solution_type tmp(x.first().size(), Zero);

            for (integer i = 0; i < A_.nCells(); ++i) {
                tmp += x[i] * y[i];
            }
            return tmp;
        }
        template <>
        void GMRES<realSquareMatrix>::update(const Array<Array<solution_type>> &H,
                                             Array<solution_type> &x, Array<solution_type> &zets,
                                             const Array<Array<solution_type>> &v,
                                             integer iDim) const {
            for (integer i = iDim; i >= 0; --i) {
                zets[i] /= H[i][i];
                for (integer j = i - 1; j >= 0; --j) {
                    zets[j] -= zets[i] * H[j][i];
                }
            }
            for (integer i = 0; i <= iDim; ++i) {
                for (integer n = 0; n < A_.nCells(); ++n) {
                    x[n] += v[i][n] * zets[i];
                }
            }
        }

        template <>
        void GMRES<realSquareMatrix>::solve(Array<solution_type> &x, const Array<solution_type> &b,
                                            const integer nIter) const {
            if (smootherPtr_) {
                smootherPtr_->solve(x, b, nFirstSolve_);
            }
            Array<solution_type> residual(A_.nCells());
            Array<solution_type> Ax(A_.nCells());
            for (auto &e : residual) {
                e.resize(x.first().size());
            }
            for (auto &e : Ax) {
                e.resize(x.first().size());
            }

            auto precond = sparseSolver::preconditioner<realSquareMatrix>::creator(A_, gmCont_);

            A_.multiple(x, Ax);
            precond->precondition(residual, b - Ax);

            auto beta = L2Norm(residual);

            if (max(beta) < epsilon_) {
                return;
            }
            for (auto &e : beta) {
                if (e == 0) {
                    return;
                }
            }

            // Allocate the space for the givens rotations, and the upper
            // Hessenburg matrix.
            Array<Array<solution_type>> H(krylovDimension_ + 1);
            for (auto &e : H) {
                e.resize(krylovDimension_);
            }

            // The Givens rotations include the sine and cosine term. The
            // cosine term is in column zero, and the sine term is in column
            // one.
            Array<Array<solution_type>> givens(krylovDimension_ + 1);
            for (auto &e : givens) {
                e.resize(2);
                for (auto &ee : e) {
                    ee.resize(x.first().size());
                }
            }

            Array<Array<solution_type>> v(krylovDimension_ + 1);
            for (auto &e : v) {
                e.resize(A_.nCells());
                for (auto &ee : e) {
                    ee.resize(x.first().size());
                }
            }

            Array<solution_type> zeta(krylovDimension_ + 1);
            for (auto &ee : zeta) {
                ee.resize(x.first().size());
            }

            solution_type relres(x.first().size(), Zero);
            for (integer k = 0; k < nRestarts_; ++k) {

                for (integer ii = 0; ii < residual.size(); ++ii) {
                    v[0][ii] = residual[ii] / beta;
                }
                zeta.setZero();
                zeta[0] = beta;
                integer m = 0;
                for (integer j = 0; j < krylovDimension_; ++j) {
                    A_.multiple(v[j], Ax);
                    precond->precondition(v[j + 1], Ax);

                    m = j + 1;
                    for (integer i = 0; i <= j; ++i) {
                        H[i][j] = dotPro(v[i], v[j + 1]);
                        for (integer ii = 0; ii < residual.size(); ++ii) {
                            v[j + 1][ii] -= (H[i][j] * v[i][ii]);
                        }
                    }

                    H[j + 1][j] = L2Norm(v[j + 1]);

                    for (auto &e : H[j + 1][j]) {
                        if (e == 0) {
                            break;
                        }
                    }

                    for (integer ii = 0; ii < residual.size(); ++ii) {
                        v[j + 1][ii] /= H[j + 1][j];
                    }

                    for (integer i = 0; i < j; ++i) {
                        auto tmp = givens[i][0] * H[i][j] + givens[i][1] * H[i + 1][j];
                        H[i + 1][j] = -givens[i][1] * H[i][j] + givens[i][0] * H[i + 1][j];
                        H[i][j] = tmp;
                    }

                    for (integer nn = 0; nn < x.first().size(); ++nn) {
                        if (H[j + 1][j][nn] == 0) {
                            givens[j][0][nn] = 1;
                            givens[j][1][nn] = 0;
                        } else if (mag(H[j][j][nn]) > mag(H[j + 1][j][nn])) {
                            auto tau = H[j + 1][j][nn] / H[j][j][nn];
                            givens[j][0][nn] = inv(sqrt(real(1) + sqr(tau)));
                            givens[j][1][nn] = givens[j][0][nn] * tau;
                        } else {
                            auto tau = H[j][j][nn] / H[j + 1][j][nn];
                            givens[j][1][nn] = inv(sqrt(real(1) + sqr(tau)));
                            givens[j][0][nn] = givens[j][1][nn] * tau;
                        }
                    }

                    solution_type hjj = H[j][j];
                    solution_type hjp1j = H[j + 1][j];

                    H[j][j] = givens[j][0] * hjj + givens[j][1] * hjp1j;
                    H[j + 1][j] = -givens[j][1] * hjj + givens[j][0] * hjp1j;

                    solution_type zj = zeta[j];
                    solution_type zjp1 = zeta[j + 1];
                    zeta[j] = givens[j][0] * zj + givens[j][1] * zjp1;
                    zeta[j + 1] = -givens[j][1] * zj + givens[j][0] * zjp1;

                    relres = mag(zeta[j + 1]);
                    bool allConverg = true;
                    for (integer nn = 0; nn < x.first().size(); ++nn) {
                        if (relres[nn] >= epsilon_ * beta[nn]) {
                            allConverg = false;
                            break;
                        }
                    }
                    if (allConverg) {
                        break;
                    }
                }

                update(H, x, zeta, v, m - 1);

                A_.multiple(x, Ax);
                precond->precondition(residual, b - Ax);

                beta = L2Norm(residual);
                for (auto &e : beta) {
                    if (e == 0) {
                        return;
                    }
                }
            }
        }
    } // namespace sparseSolver
} // namespace OpenHurricane