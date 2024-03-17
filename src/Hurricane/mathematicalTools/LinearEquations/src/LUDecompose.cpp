/*!
 * \file LUDecompose.cpp
 * \brief Main subroutines for the LU decompose.
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

#include "LUDecompose.hpp"
#include "Lists.hpp"

// Doolittle LU decompose
void OpenHurricane::LUDecomposeD(realSquareMatrix &A) {
    integer n = A.n();

    for (integer k = 0; k < n; ++k) {
        for (integer j = k; j < n; ++j) {
            real sumlu = Zero;
            for (integer t = 0; t < k; ++t) {
                sumlu += A(k, t) * A(t, j);
            }
            A(k, j) -= sumlu;
        }

        for (integer i = k + 1; i < n; ++i) {
            if (k < n - 1) {
                real sumlu = Zero;
                for (integer t = 0; t < k; ++t) {
                    sumlu += A(i, t) * A(t, k);
                }
                A(i, k) = (A(i, k) - sumlu) / A(k, k);
            }
        }
    }
}

// Doolittle LU decompose with pivoting
void OpenHurricane::LUDecomposeDPivoting(realSquareMatrix &A, integerList &M) {
    integer sign;
    LUDecomposeDPivoting(A, M, sign);
}

// Crout LU decompose
void OpenHurricane::LUDecomposeC(realSquareMatrix &A) {
    integer n = A.n();

    for (integer k = 0; k < n; ++k) {
        for (integer i = k; i < n; ++i) {
            real sumlu = Zero;
            for (integer t = 0; t < k; ++t) {
                sumlu += A(i, t) * A(t, k);
            }
            A(i, k) -= sumlu;
        }

        for (integer j = k + 1; j < n; ++j) {
            if (k < n - 1) {
                real sumlu = Zero;
                for (integer t = 0; t < k; ++t) {
                    sumlu += A(k, t) * A(t, j);
                }
                A(k, j) = (A(k, j) - sumlu) / A(k, k);
            }
        }
    }
}

void OpenHurricane::LUDecomposeDBanded(realSquareMatrix &A, const integer s,
                                       const integer r) {
    integer n = A.n();

    for (integer k = 0; k < n; ++k) {
        for (integer j = k; j < min(k + s + 1, n); ++j) {
            real sumlu = Zero;
            for (integer t = max(0, max(k - r, k - r)); t < k; ++t) {
                sumlu += A(k, t) * A(t, j);
            }
            A(k, j) -= sumlu;
        }

        for (integer i = k + 1; i < min(k + r + 1, n); ++i) {
            if (k < n - 1) {
                real sumlu = Zero;
                for (integer t = max(i - r, max(k - s, 0)); t < k; ++t) {
                    sumlu += A(i, t) * A(t, k);
                }
                A(i, k) = (A(i, k) - sumlu) / A(k, k);
            }
        }
    }
}
