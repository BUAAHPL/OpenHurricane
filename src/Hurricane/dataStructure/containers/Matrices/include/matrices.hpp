/*!
 * \file matrices.hpp
 * \brief Header of matrices
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
#include "DiagonalMatrix.hpp"
#include "RectangularMatrix.hpp"
#include "SquareMatrix.hpp"
#include "SymmMatrix.hpp"
#include "realArray.hpp"

namespace OpenHurricane {
    using realSquareMatrix = SquareMatrix<real>;
    using realRectangularMatrix = RectangularMatrix<real>;
    using realDiagonalMatrix = DiagonalMatrix<real>;
    using realSymmMatrix = SymmMatrix<real>;

    /*!\brief Return the inverse matrix.*/
    hur_nodiscard realSquareMatrix inv(const realSquareMatrix &M);

    /*!\brief Return the inverse matrix.*/
    hur_nodiscard realSquareMatrix inv(const realSymmMatrix &M);

    hur_nodiscard inline realSquareMatrix operator/(const realSquareMatrix &M,
                                                    const realSquareMatrix &N) {
        return inv(N) * M;
    }

    hur_nodiscard inline realSquareMatrix toInverse(const realSquareMatrix &M) {
        return inv(M);
    }

    template <> hur_nodiscard SquareMatrix<real> SquareMatrix<real>::inv() const;

    /*!\brief Return the the eigenvalues of real symmetry matrix.*/
    hur_nodiscard realArray eigenValues(const realSymmMatrix &M);

    /**
     * \brief Transpose the square matrix itself.
     */
    void transposing(realSquareMatrix &M);

    /**
     * \brief Transpose the square matrix itself.
     */
    inline void transposing(realSymmMatrix &M) {
        // Do nothing
    }

} //  namespace OpenHurricane
