/*!
 * \file pointArrays.hpp
 * \brief Headers of the geometry array based on point.
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

#include "geometryArrays.hpp"
#include "pointMesh.hpp"

namespace OpenHurricane {
    using pointRealArray = realGeometryArray<pointMesh>;
    using pointComplexArray = complexGeometryArray<pointMesh>;
    using pointVector2DArray = vector2DGeometryArray<pointMesh>;
    using pointVectorArray = vectorGeometryArray<pointMesh>;
    using pointTensorArray = tensorGeometryArray<pointMesh>;
    using pointSymmTensorArray = symmTensorGeometryArray<pointMesh>;

    template <>
    hur_nodiscard inline realArray
    geometryArray<real, pointMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    template <>
    hur_nodiscard inline realArray
    geometryArray<vector2D, pointMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    template <>
    hur_nodiscard inline realArray
    geometryArray<vector, pointMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    template <>
    hur_nodiscard inline realArray
    geometryArray<tensor, pointMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    template <>
    hur_nodiscard inline realArray
    geometryArray<symmTensor, pointMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
} // namespace OpenHurricane
