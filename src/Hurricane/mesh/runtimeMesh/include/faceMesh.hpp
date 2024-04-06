/*!
 * \file faceMesh.hpp
 * \brief Headers of the face mesh.
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

#include "boundary.hpp"
#include "runtimeMesh.hpp"
namespace OpenHurricane {
    /*!\brief The class of face mesh.*/
    class faceMesh {
    public:
        using Mesh = runtimeMesh;

    private:
        const Mesh &mesh_;

    public:
        hur_nodiscard static const OpenHurricane::faceMesh &nullObject() {
            return NullRefObj::nullRef<faceMesh>();
        }

        faceMesh() : mesh_(Mesh::nullObject()) {}

        explicit faceMesh(const Mesh &mesh) : mesh_(mesh) {}

        hur_nodiscard inline integer size() const noexcept { return mesh_.nTotalFaces(); }

        static hur_nodiscard integer size(const Mesh &mesh) noexcept { return mesh.nTotalFaces(); }

        hur_nodiscard inline integer internalArraySize() const noexcept { return mesh_.nFaces(); }

        static hur_nodiscard integer internalArraySize(const Mesh &mesh) noexcept {
            return mesh.nFaces();
        }

        hur_nodiscard const Mesh &operator()() const noexcept { return mesh_; }
    };
} // namespace OpenHurricane