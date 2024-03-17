/*!
 * \file limitersForLinear.hpp
 * \brief Header of limiter for piecewise linear reconstruction.
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
#include "cellArrays.hpp"
#include "controller.hpp"
#include "faceArrays.hpp"
#include "objectFactory.hpp"
#include "smartPointer.hpp"

namespace OpenHurricane {
    /*!\brief The base class of limiters for piecewise linear reconstruction.*/
    class limitersForLinear {
    public:
        declareClassNames;

        declareObjFty(limitersForLinear, controller, (const controller &cont), (cont));

        inline limitersForLinear() {}

        /*!\brief Select null constructed.*/
        static uniquePtr<limitersForLinear> creator(const controller &cont);

        /*!\brief Destructor.*/
        virtual ~limitersForLinear() noexcept {}

        /*!\brief Calculating the limiters.*/
        virtual void calcLimiters(
            const geometryArray<real, cellMesh> &cellQ,
            const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &gradQ,
            geometryArray<real, cellMesh> &_limiters) const = 0;

        /*!\brief Calculating the limiters.*/
        virtual void calcLimiters(
            const geometryArray<vector, cellMesh> &cellQ,
            const geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &gradQ,
            geometryArray<vector, cellMesh> &_limiters) const = 0;

        void updateBoundaryField(geometryArray<real, cellMesh> &_limiters) const;
        void updateBoundaryField(geometryArray<vector, cellMesh> &_limiters) const;

        /*!\brief Calculating the limiters and updating the boundary field.*/
        template <class Type>
        void
        limiters(const geometryArray<Type, cellMesh> &cellQ,
                 const geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &gradQ,
                 geometryArray<Type, cellMesh> &_limiters) const;

        /*!\brief Calculating the limiters and updating the boundary field.*/
        template <class Type> void limiters(geometryArray<Type, cellMesh> &cellQ) const;

        /*!\brief Calculating the limiters and updating the boundary field.*/
        template <class Type>
        geometryArray<Type, cellMesh> limiters(
            const geometryArray<Type, cellMesh> &cellQ,
            const geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &gradQ) const;
    };

} // namespace OpenHurricane

#include "limitersForLinear.inl"
