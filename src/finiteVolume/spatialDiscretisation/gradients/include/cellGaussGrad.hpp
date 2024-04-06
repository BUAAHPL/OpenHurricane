/*!
 * \file cellGaussGrad.hpp
 * \brief Headers of the cell-based Green-Gauss gradient approach.
 *        The subroutines and functions are in the <i>cellGaussGrad.cpp</i> file.
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
#include "gradient.hpp"

namespace OpenHurricane {

    class cellGaussGrad : public gradient {
    public:
        declareClassNames;

        cellGaussGrad(const controller &cont) {}

        virtual inline ~cellGaussGrad() noexcept {}

        void
        gaussGrad(const geometryArray<real, cellMesh> &cellQ,
                  geometryArray<typename outerProduct<vector, real>::type, cellMesh> &grad) const;

        void
        gaussGrad(const geometryArray<vector, cellMesh> &cellQ,
                  geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &grad) const;

        virtual void calcGrad(const geometryArray<real, cellMesh> &cellQ,
                              geometryArray<vector, cellMesh> &grad) const;

        virtual void calcGrad(const geometryArray<vector, cellMesh> &cellQ,
                              geometryArray<tensor, cellMesh> &grad) const;
    };

} // namespace OpenHurricane

inline void OpenHurricane::cellGaussGrad::calcGrad(const geometryArray<real, cellMesh> &cellQ,
                                                   geometryArray<vector, cellMesh> &grad) const {
    gaussGrad(cellQ, grad);
}

inline void OpenHurricane::cellGaussGrad::calcGrad(const geometryArray<vector, cellMesh> &cellQ,
                                                   geometryArray<tensor, cellMesh> &grad) const {
    gaussGrad(cellQ, grad);
}
