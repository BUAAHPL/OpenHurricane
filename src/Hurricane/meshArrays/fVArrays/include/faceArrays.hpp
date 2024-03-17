/*!
 * \file faceArrays.hpp
 * \brief Headers of the geometry field based on face.
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

#include "faceMesh.hpp"
#include "geometryArrays.hpp"

namespace OpenHurricane {
    using faceRealArray = realGeometryArray<faceMesh>;
    using faceComplexArray = complexGeometryArray<faceMesh>;
    using faceVector2DArray = vector2DGeometryArray<faceMesh>;
    using faceVectorArray = vectorGeometryArray<faceMesh>;
    using faceTensorArray = tensorGeometryArray<faceMesh>;
    using faceSymmTensorArray = symmTensorGeometryArray<faceMesh>;

    template <>
    void geometryArray<real, faceMesh>::writeOutput(fileOsstream &fos, const integer fzid) const;

    template <>
    void geometryArray<vector, faceMesh>::writeOutput(fileOsstream &fos, const integer fzid) const;

    template <>
    void geometryArray<vector2D, faceMesh>::writeOutput(fileOsstream &fos,
                                                        const integer fzid) const;

    template <>
    void geometryArray<tensor, faceMesh>::writeOutput(fileOsstream &fos, const integer fzid) const;

    template <>
    void geometryArray<symmTensor, faceMesh>::writeOutput(fileOsstream &fos,
                                                          const integer fzid) const;

    template <>
    void geometryArray<real, faceMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                          const integer fzid) const;

    template <>
    void geometryArray<vector, faceMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                            const integer fzid) const;

    template <>
    void geometryArray<vector2D, faceMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                              const integer fzid) const;

    template <>
    void geometryArray<tensor, faceMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                            const integer fzid) const;

    template <>
    void geometryArray<symmTensor, faceMesh>::writeMinMaxOutput(fileOsstream &fos,
                                                                const integer fzid) const;
    template <>
    hur_nodiscard inline realArray geometryArray<real, faceMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    template <>
    hur_nodiscard inline realArray
    geometryArray<vector2D, faceMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    template <>
    hur_nodiscard inline realArray
    geometryArray<vector, faceMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    template <>
    hur_nodiscard inline realArray
    geometryArray<tensor, faceMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
    template <>
    hur_nodiscard inline realArray
    geometryArray<symmTensor, faceMesh>::realComponent(const int i) const {
        return Base::component(i);
    }
} // namespace OpenHurricane