/*!
 * \file gradient.hpp
 * \brief Headers of the gradient on face.
 *        The subroutines and functions are in the <i>faceGrad.cpp</i> file.
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
#include "faceArrays.hpp"
#include "runtimeMesh.hpp"

namespace OpenHurricane {
    /*!\brief The class of face gradient.*/
    class faceGrad {
    private:
    public:
        /*!\brief Construct from components.*/
        faceGrad();

        static typename outerProduct<vector, real>::type
        grad(const geometryArray<real, cellMesh> &cellQ,
             const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &cellGrad,
             const integer faceI);

        static void
        grad(const geometryArray<real, cellMesh> &cellQ,
             const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &cellGrad,
             geometryArray<typename outerProduct<vector, real>::type, faceMesh> &faceGrads);

        static typename outerProduct<vector, vector>::type
        grad(const geometryArray<vector, cellMesh> &cellQ,
             const geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &cellGrad,
             const integer faceI);

        static void
        grad(const geometryArray<vector, cellMesh> &cellQ,
             const geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &cellGrad,
             geometryArray<typename outerProduct<vector, vector>::type, faceMesh> &faceGrads);

        static geometryArray<typename outerProduct<vector, real>::type, faceMesh>
        grad(const geometryArray<real, cellMesh> &cellQ,
             const geometryArray<typename outerProduct<vector, real>::type, cellMesh> &cellGrad);
        static geometryArray<typename outerProduct<vector, vector>::type, faceMesh>
        grad(const geometryArray<vector, cellMesh> &cellQ,
             const geometryArray<typename outerProduct<vector, vector>::type, cellMesh> &cellGrad);
    };

} // namespace OpenHurricane

#include "faceGrad.inl"