/*!
 * \file leastSquareGrad.hpp
 * \brief Headers of the least square gradient approach.
 *        The subroutines and functions are in the <i>leastSquareGrad.cpp</i> file.
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

    class leastSquareGrad : public gradient {
    private:
        string key_;

    public:
        declareClassNames;

        leastSquareGrad(const controller &cont)
            : key_(cont.findWordOrDefault("weightType", "WLSQ2")) {}

        inline virtual ~leastSquareGrad() noexcept {}

        virtual void calcGrad(const geometryArray<real, cellMesh> &cellQ,
                              geometryArray<vector, cellMesh> &grad) const;

        virtual void calcGrad(const geometryArray<vector, cellMesh> &cellQ,
                              geometryArray<tensor, cellMesh> &grad) const;

        template <class Type>
        void
        WLSGrad(const geometryArray<Type, cellMesh> &cellQ,
                geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &grad) const;
    };

} // namespace OpenHurricane
