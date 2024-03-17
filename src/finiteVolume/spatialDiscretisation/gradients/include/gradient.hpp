/*!
 * \file gradient.hpp
 * \brief Headers of the gradient.
 *        The subroutines and functions are in the <i>gradient.cpp</i> file.
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
#include "runtimeMesh.hpp"

namespace OpenHurricane {

    /*!\brief The base class of gradient scheme.*/
    class gradient {
    public:
        // Static data

        declareClassNames;
        declareObjFty(gradient, controller, (const controller &cont), (cont));

        static uniquePtr<gradient> creator(const controller &cont);

        inline virtual ~gradient() noexcept {}

        virtual void calcGrad(const geometryArray<real, cellMesh> &cellQ,
                              geometryArray<vector, cellMesh> &grad) const = 0;

        virtual void calcGrad(const geometryArray<vector, cellMesh> &cellQ,
                              geometryArray<tensor, cellMesh> &grad) const = 0;

        template <class Type>
        void grad(const geometryArray<Type, cellMesh> &cellQ,
                  geometryArray<typename outerProduct<vector, Type>::type, cellMesh> &grad) const;

        template <class Type>
        OpenHurricane::geometryArray<typename outerProduct<vector, Type>::type, cellMesh>
        grad(const geometryArray<Type, cellMesh> &cellQ, const string &name) const;

        template <class Type> void grad(geometryArray<Type, cellMesh> &cellQ) const;

        static void updateBoundaryField(geometryArray<vector, cellMesh> &grad);
        static void updateBoundaryField(geometryArray<tensor, cellMesh> &grad);
    };

} // namespace OpenHurricane

#include "gradient.inl"