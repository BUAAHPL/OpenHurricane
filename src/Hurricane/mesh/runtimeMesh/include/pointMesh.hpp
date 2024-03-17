/*!
 * \file pointMesh.hpp
 * \brief Headers of the point mesh.
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

#include "runtimeMesh.hpp"

namespace OpenHurricane {
    /*!\brief The class of point mesh.*/
    class pointMesh {
    public:
        using Mesh = runtimeMesh;

    private:
        const Mesh &mesh_;

    public:
        hur_nodiscard static const OpenHurricane::pointMesh &nullObject() {
            return NullRefObj::nullRef<pointMesh>();
        }

        pointMesh() : mesh_(Mesh::nullObject()) {}

        explicit pointMesh(const Mesh &mesh) : mesh_(mesh) {}

        hur_nodiscard inline integer size() const noexcept { return mesh_.nTotalPoints(); }

        hur_nodiscard static integer size(const Mesh &mesh) noexcept { return mesh.nTotalPoints(); }

        hur_nodiscard inline integer internalArraySize() const noexcept { return mesh_.nPoints(); }

        hur_nodiscard static integer internalArraySize(const Mesh &mesh) noexcept {
            return mesh.nPoints();
        }

        hur_nodiscard inline const Mesh &operator()() const noexcept { return mesh_; }
    };
} // namespace OpenHurricane