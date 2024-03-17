/*!
 * \file meshInterpolation.hpp
 * \brief Header of meshInterpolation.
 *       The subroutines and functions are in the <i>meshInterpolation.cpp</i> file.
 * \author Yang Hongzhen
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

#include "ArrayInclude.hpp"
#include "OpenHurricane.hpp"
#include "LUDecompose.hpp"
#include "matrices.hpp"

namespace OpenHurricane {

    class meshInterpolation {
    private:
        /*!\brief Interpolation base points.*/
        const integerArrayArray &nbr_;

        /*!\brief Coordinate interpolation source points.*/
        const vectorArray &sNodes_;

        /*!\brief Coordinate interpolation target points.*/
        const vectorArray &tNodes_;

        ///*!\brief Interploation matrix.*/
        Array<realSymmMatrix> M_;

        integerArrayArray pivotIndices_;

    public:
        meshInterpolation(const integerArrayArray &nbr, const vectorArray &sNode,
                          const vectorArray &tNode);

        template <class Type> void interpolate(const Array<Type> &sor, Array<Type> &tar);

    private:
        void setM();

        hur_nodiscard inline real fai(const real x) const noexcept { return 1.0 / (1.0 + x * x); }
    };

    template <> void meshInterpolation::interpolate(const Array<real> &sor, Array<real> &tar);
} // namespace OpenHurricane

template <class Type>
void OpenHurricane::meshInterpolation::interpolate(const Array<Type> &sor, Array<Type> &tar) {
    setM();
    tar = Zero;
    for (integer n = 0; n < nbr_.size(); n++) {
        const integer k = nbr_[n].size();
        for (integer i = 0; i < feature<Type>::nElements_; i++) {
            realArray D(k, Zero);
            for (integer j = 0; j < k; j++) {
                const integer pi = nbr_[n][j];
                D[j] = sor[pi][i];
            }
            LUBacksubstituteDPivoting(M_[n], D, pivotIndices_[n]);
            for (integer j = 0; j < k; j++) {
                const integer pi = nbr_[n][j];
                const real rd = (tNodes_[n] - sNodes_[pi]).magnitude();
                tar[n][i] += D[j] * fai(rd);
            }
        }
    }
}