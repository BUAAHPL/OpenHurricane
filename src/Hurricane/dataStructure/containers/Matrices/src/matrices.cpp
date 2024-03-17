/*!
 * \file matrices.cpp
 * \brief Main subroutines of matrices.
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
#include "matrices.hpp"
#include "LUDecompose.hpp"

OpenHurricane::realSquareMatrix OpenHurricane::inv(const realSquareMatrix &A) {
    if (A.onlyDiagFillNonZero()) {
        realSquareMatrix invA(A);
        for (realSquareMatrix::size_type i = 0; i < A.n(); ++i) {
            invA(i, i) = inv(A(i, i));
        }
        invA.setOnlyDiagFillNonZero();
        return invA;
    } else {
        realSquareMatrix invA(A.m(), A.n(), Zero);
        realSquareMatrix AA = A;
        integerList M(A.m());
        LUDecomposeDPivoting(AA, M);
        List<real> bb(A.n(), Zero);
        for (integer i = 0; i < A.n(); ++i) {
            bb[i] = 1.0;
            LUBacksubstituteDPivoting(AA, bb, M);
            for (integer j = 0; j < A.m(); ++j) {
                invA(i, j) = bb[j];
            }
            bb = Zero;
        }
        transposing(invA);
        return invA;
    }
}

OpenHurricane::realSquareMatrix OpenHurricane::inv(const realSymmMatrix &M) {
    return inv(dynamic_cast<const realSquareMatrix &>(M));
}

template <>
hur_nodiscard OpenHurricane::SquareMatrix<OpenHurricane::real>
OpenHurricane::SquareMatrix<OpenHurricane::real>::inv() const {
    return toInverse(*this);
}

OpenHurricane::realArray OpenHurricane::eigenValues(const realSymmMatrix &M) {
    integer n = M.m();
    integer nRow = 1;
    integer nCol = 0;
    real apq = M[0][1];
    realSymmMatrix A(M);
    realArray eigenvalues(n, Zero);
    integer counter = 0;
    real precision = 1.e-10;
    integer maxStep = 1000000;
    while (true) {
        // Step 1: Found the largest element on the off-diagnoal line of the matrix M.
        for (integer i = 0; i < n; i++) {
            for (integer j = 0; j < n; j++) {
                if (i == j)
                    continue;
                else {
                    real d = fabs(A[i][j]);
                    if (d > apq) {
                        apq = d;
                        nRow = i;
                        nCol = j;
                    }
                }
            }
        }
        counter++;
        // Step 2: Percision check
        if (apq < precision)
            break;
        else {
            if (counter > maxStep) {
                LFatal("Cannot gain effect eigenvalue solution under current precision limits");
            }
        }

        // Step 3: Calculate rotate angle
        real app = A[nRow][nRow];
        real aqq = A[nCol][nCol];
        real psi = 0.5 * atan2(-2 * apq, aqq - app);
        real cospsi = cos(psi);
        real sinpsi = sin(psi);
        real cos2psi = cos(2 * psi);
        real sin2psi = sin(2 * psi);

        // Step 4: Update new matrix
        A[nRow][nRow] = app * cospsi * cospsi + aqq * sinpsi * sinpsi + 2 * apq * cospsi * sinpsi;
        A[nCol][nCol] = app * sinpsi * sinpsi + aqq * cospsi * cospsi - 2 * apq * cospsi * sinpsi;
        A[nRow][nCol] = 0.5 * (aqq - app) * sin2psi + apq * cos2psi;
        A[nCol][nRow] = A[nRow][nCol];
        for (integer i = 0; i < n; i++) {
            if (i == nRow || i == nCol) {
                continue;
            } else {
                A[nRow][i] = cospsi * A[nRow][i] + sinpsi * A[nCol][i];
                A[nCol][i] = cospsi * A[nCol][i] - sinpsi * A[nRow][i];
            }
        }
        for (integer j = 0; j < n; j++) {
            if (j == nRow || j == nCol) {
                continue;
            } else {
                A[j][nRow] = cospsi * A[j][nRow] + sinpsi * A[j][nCol];
                A[j][nCol] = cospsi * A[j][nCol] - sinpsi * A[j][nRow];
            }
        }
        for (integer k = 0; k < n; k++) {
            eigenvalues[k] = A[k][k];
        }
    }
    return eigenvalues;
}

void OpenHurricane::transposing(realSquareMatrix &M) {
    for (integer i = 0; i < M.m(); ++i) {
        for (integer j = i + 1; j < M.n(); ++j) {
            real tmp = M(i, j);
            M(i, j) = M(j, i);
            M(j, i) = tmp;
        }
    }
}
