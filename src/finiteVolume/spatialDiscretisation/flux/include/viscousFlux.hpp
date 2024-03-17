/*!
 * \file viscousFlux.hpp
 * \brief Headers of viscous flux.
 *        The subroutines and functions are in the <i>viscousFlux.cpp</i> file.
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
#include "faceGrad.hpp"

#include "faceInterpolation.hpp"
#include "gradients.hpp"

namespace OpenHurricane {
    template <class Type>
    static void
    visFluxPerZone(const geometryArray<typename outerProduct<vector, Type>::type, faceMesh> &gf,
                   const integer faceZoneId, Array<Type> &rhs);

    void visFluxPerZone(const geometryArray<symmTensor, faceMesh> &gf, const integer faceZoneId,
                        Array<vector> &rhs);

    template <class Type>
    static Array<Type>
    visFlux(const geometryArray<typename outerProduct<vector, Type>::type, faceMesh> &gf);

    template <class Type>
    static void
    visFlux(const geometryArray<typename outerProduct<vector, Type>::type, faceMesh> &gf,
            geometryArray<Type, cellMesh> &cellQ);

    Array<vector> visFlux(const geometryArray<symmTensor, faceMesh> &gf);

    void visFlux(const geometryArray<symmTensor, faceMesh> &gf,
                 geometryArray<vector, cellMesh> &cellQ);
} // namespace OpenHurricane

#include "viscousFlux.inl"
