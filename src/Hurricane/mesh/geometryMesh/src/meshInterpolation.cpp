/*!
 * \file meshInterpolation.cpp
 * \brief The subroutines and functions of CFD time advance iteration
 * \author Yang Hongzhen
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

#include "meshInterpolation.hpp"

OpenHurricane::meshInterpolation::meshInterpolation(const integerArrayArray &nbr,
                                                    const vectorArray &sNode,
                                                    const vectorArray &tNode)
    : nbr_(nbr), sNodes_(sNode), tNodes_(tNode), M_(), pivotIndices_() {}

void OpenHurricane::meshInterpolation::setM() {
    if (M_.size() == 0) {
        M_.resize(nbr_.size());
        pivotIndices_.resize(nbr_.size());

        for (integer n = 0; n < nbr_.size(); n++) {
            const integer k = nbr_[n].size();
            M_[n].resize(k);
            pivotIndices_[n].resize(k);

            for (integer i = 0; i < k; i++) {
                const integer pi = nbr_[n][i];
                for (integer j = i; j < k; j++) {
                    const integer pj = nbr_[n][j];
                    const real rd = (sNodes_[pi] - sNodes_[pj]).magnitude();
                    M_[n][i][j] = fai(rd);
                }
            }
            LUDecomposeDPivoting(M_[n], pivotIndices_[n]);
        }
    }
}

template <>
void OpenHurricane::meshInterpolation::interpolate(const Array<real> &sor,
                                                   Array<real> &tar) {
    setM();
    tar = Zero;
    for (integer n = 0; n < nbr_.size(); n++) {
        const integer k = nbr_[n].size();
        realArray D(k, Zero);
        for (integer i = 0; i < k; i++) {
            const integer pi = nbr_[n][i];
            D[i] = sor[pi];
        }

        LUBacksubstituteDPivoting(M_[n], D, pivotIndices_[n]);
        for (integer i = 0; i < k; i++) {
            const integer pi = nbr_[n][i];
            const real rd = (tNodes_[n] - sNodes_[pi]).magnitude();
            tar[n] += D[i] * fai(rd);
        }
    }
}
