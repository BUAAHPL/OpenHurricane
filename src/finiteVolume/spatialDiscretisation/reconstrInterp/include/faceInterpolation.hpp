/*!
 * \file faceInterpolation.hpp
 * \brief Headers of the face interpolation.
 *        The subroutines and functions are in the <i>faceInterpolation.cpp</i> file.
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

#include "cellArrays.hpp"
#include "faceArrays.hpp"
#include "gradients.hpp"
#include "transfer.hpp"

namespace OpenHurricane {
    namespace fv {
        /*!\brief Interpolation field onto face.*/
        template <class Type>
        geometryArray<Type, faceMesh> interpolate(const geometryArray<Type, cellMesh> &cf);

        /*!\brief Interpolation field onto face.*/
        template <>
        geometryArray<real, faceMesh> interpolate(const geometryArray<real, cellMesh> &cf);

        /*!\brief Interpolation field onto face.*/
        template <>
        geometryArray<vector, faceMesh> interpolate(const geometryArray<vector, cellMesh> &cf);

        /*!\brief Interpolation field onto face.*/
        template <class Type>
        void interpolate(const geometryArray<Type, cellMesh> &cf,
                         geometryArray<Type, faceMesh> &ff);

        /*!\brief Interpolation field onto face zone by given zone id.*/
        template <class Type>
        Array<Type> interpolate(const geometryArray<Type, cellMesh> &cf, const integer faceZoneId);

        /*!\brief Interpolation field onto face zone by given zone id.*/
        template <class Type>
        Array<Type> interpolateNoCorr(const geometryArray<Type, cellMesh> &cf,
                                      const integer faceZoneId);

        /**\brief Interpolation field onto ghost cell by given boundary face value.*/
        void extrapolate(const Array<Array<real>> &faceVar, geometryArray<real, cellMesh> &var);

        /**\brief Interpolation field onto ghost cell by given boundary face value.*/
        void extrapolate(const Array<Array<vector>> &faceVar, geometryArray<vector, cellMesh> &var);

        /**\brief Interpolation field onto ghost cell by given boundary face value and face zone index.*/
        template <class Type>
        void extrapolate(const Array<Type> &faceVar, geometryArray<Type, cellMesh> &var,
                         const integer faceZoneId);
    } // namespace fv
} // namespace OpenHurricane

#include "faceInterpolation.inl"