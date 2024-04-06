/*!
 * \file GaussElimination.hpp
 * \brief Headers of the Gauss elimination.
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

#include "matrices.hpp"

namespace OpenHurricane {
    // Gauss elimination (will change A and b)
    template <class Type> void GaussElimination(realSquareMatrix &A, List<Type> &b) {
#ifdef HUR_DEBUG
        if (A.m() != b.size()) {
            LFatal("A.m!=b.size()");
        }
#endif // HUR_DEBUG

        integer n = A.m();

        // Elimination
        for (integer k = 0; k < n - 1; ++k) {
            for (integer i = k + 1; i < n; ++i) {
                if (mag(A(k, k)) < 1e-20) {
                    LFatal("Gauss elimination failed");
                }
                Type Mik = A(i, k) / A(k, k);
                for (integer j = k + 1; j < n; ++j) {
                    //A(i, j) = A(i, j) - Mik * A(k, j);
                    A(i, j) -= Mik * A(k, j);
                }
                //b[i] = b[i] - Mik * b[k];
                b[i] -= Mik * b[k];
            }
        }

        // Back-substitution
        b[n - 1] = b[n - 1] / A(n - 1, n - 1);

        for (integer k = n - 2; k >= 0; --k) {
            Type sumax = Zero;
            for (integer j = k + 1; j < n; ++j) {
                sumax += A(k, j) * b[j];
            }
            b[k] = (b[k] - sumax) / A(k, k);
        }
    }

    // Gauss elimination with pivoting (will change A and b)
    template <class Type> void GaussEliminationPivoting(realSquareMatrix &A, List<Type> &b) {
#ifdef HUR_DEBUG
        if (A.m() != b.size()) {
            LFatal("A.m!=b.size()");
        }
#endif // HUR_DEBUG

        integer n = A.m();

        // Elimination
        for (integer k = 0; k < n - 1; ++k) {
            integer ikmax = k;
            real largestA = mag(A(k, k));
            for (integer i = k + 1; i < n; ++i) {
                if (mag(A(i, k)) > largestA) {
                    ikmax = i;
                    largestA = mag(A(i, k));
                }
            }
            if (k != ikmax) {
                for (integer j = k; j < n; ++j) {
                    Swap(A(k, j), A(ikmax, j));
                }
                Swap(b[k], b[ikmax]);
            }
            for (integer i = k + 1; i < n; ++i) {
                if (mag(A(k, k)) < 1e-20) {
                    LFatal("Gauss elimination failed");
                }
                Type Mik = A(i, k) / A(k, k);
                for (integer j = k + 1; j < n; ++j) {
                    //A(i, j) = A(i, j) - Mik * A(k, j);
                    A(i, j) -= Mik * A(k, j);
                }
                //b[i] = b[i] - Mik * b[k];
                b[i] -= Mik * b[k];
            }
        }

        // Back-substitution
        b[n - 1] = b[n - 1] / A(n - 1, n - 1);

        for (integer k = n - 2; k >= 0; --k) {
            Type sumax = Zero;
            for (integer j = k + 1; j < n; ++j) {
                sumax += A(k, j) * b[j];
            }
            b[k] = (b[k] - sumax) / A(k, k);
        }
    }

    // Gauss elimination (will not change A and b)
    // Return x
    template <class Type>
    void GaussElimination(const realSquareMatrix &A, const List<Type> &b, List<Type> &x) {
        realSquareMatrix tmpA = A;

        if (x.size() != b.size()) {
            x.resize(b.size());
        }
        x = b;
        GaussElimination(tmpA, x);
    }

    // Gauss elimination with pivoting (will not change A and b)
    template <class Type>
    void GaussEliminationPivoting(const realSquareMatrix &A, const List<Type> &b, List<Type> &x) {
        realSquareMatrix tmpA = A;

        if (x.size() != b.size()) {
            x.resize(b.size());
        }
        x = b;
        GaussEliminationPivoting(tmpA, x);
    }
} // namespace OpenHurricane