/*!
 * \file LUDecompose.hpp
 * \brief Headers of the LU decompose.
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
    // Doolittle LU decompose
    void LUDecomposeD(realSquareMatrix &A);

    // Doolittle LU decompose with pivoting
    void LUDecomposeDPivoting(realSquareMatrix &A, integerList &M);

    // Crout LU decompose
    void LUDecomposeC(realSquareMatrix &A);

    // Doolittle LU decompose back-substitute
    template <class Type> void LUBacksubstituteD(realSquareMatrix &A, List<Type> &b) {
        integer n = A.n();
        for (integer i = 1; i < n; ++i) {
            Type sumyb = Zero;
            for (integer t = 0; t < i; ++t) {
                sumyb += A(i, t) * b[t];
            }
            b[i] -= sumyb;
        }

        b[n - 1] = b[n - 1] / A(n - 1, n - 1);

        for (integer i = n - 2; i >= 0; --i) {
            Type sumux = Zero;
            for (integer t = i + 1; t < n; ++t) {
                sumux += A(i, t) * b[t];
            }
            b[i] = (b[i] - sumux) / A(i, i);
        }
    }

    // Doolittle LU decompose back-substitute with pivoting
    template <class Type>
    void LUBacksubstituteDPivoting(realSquareMatrix &A, List<Type> &b, integerList &M) {
        integer n = A.n();
#ifdef HUR_DEBUG
        if (A.m() != b.size()) {
            LFatal("A.m!=b.size()");
        }
        if (A.m() != M.size()) {
            LFatal("A.m!=M.size()");
        }
#endif // HUR_DEBUG

        for (integer k = 0; k < n - 1; ++k) {
            integer t = M[k];
            if (t != k) {
                Swap(b[k], b[t]);
            }
        }
        for (integer i = 0; i < n; ++i) {
            Type sumly = Zero;
            for (integer t = 0; t < i; ++t) {
                sumly += A(i, t) * b[t];
            }
            b[i] -= sumly;
        }

        for (integer i = n - 1; i >= 0; --i) {
            Type sumux = Zero;
            for (integer t = i + 1; t < n; ++t) {
                sumux += A(i, t) * b[t];
            }
            b[i] = (b[i] - sumux) / A(i, i);
        }
    }

    // Crout LU decompose back-substitute
    template <class Type> void LUBacksubstituteC(realSquareMatrix &A, List<Type> &b) {
        integer n = A.n();

        for (integer i = 0; i < n; ++i) {
            Type sumly = Zero;
            for (integer t = 0; t < i; ++t) {
                sumly += A(i, t) * b[t];
            }
            b[i] = (b[i] - sumly) / A(i, i);
        }
        for (integer i = n - 1; i >= 0; --i) {
            Type sumux = Zero;
            for (integer t = i + 1; t < n; ++t) {
                sumux += A(i, t) * b[t];
            }
            b[i] -= sumux;
        }
    }
    /**
     * \brief Doolittle LU decompose for banded matrix.
     * \param[in,out] A - The matrix for LU decomposing
     * \param[in] s - The upper bandwidth
     * \param[in] r - The lower bandwidth
     */
    void LUDecomposeDBanded(realSquareMatrix &A, const integer s, const integer r);

    /**
     * \brief Doolittle LU decompose back-substitute for banded matrix.
     * \param[in,out] A - The LU matrix
     * \param[in] s - The upper bandwidth
     * \param[in] r - The lower bandwidth
     */
    template <class Type>
    void LUBacksubstituteDBanded(realSquareMatrix &A, List<Type> &b, const integer s,
                                 const integer r) {
        integer n = A.n();
        for (integer i = 1; i < n; ++i) {
            Type sumyb = Zero;
            for (integer t = max(0, i - r); t < i; ++t) {
                sumyb += A(i, t) * b[t];
            }
            b[i] -= sumyb;
        }

        b[n - 1] = b[n - 1] / A(n - 1, n - 1);

        for (integer i = n - 2; i >= 0; --i) {
            Type sumux = Zero;
            for (integer t = i + 1; t < min(i + s + 1, n); ++t) {
                sumux += A(i, t) * b[t];
            }
            b[i] = (b[i] - sumux) / A(i, i);
        }
    }
} // namespace OpenHurricane
