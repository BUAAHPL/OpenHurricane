/*!
 * \file runtimeMesh.hpp
 * \brief Headers of the runtime mesh.
 *        The subroutines and functions are in the <i>runtimeMesh.cpp</i> file.
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
#include "geometryMesh.hpp"
#include "markRegion.hpp"

namespace OpenHurricane {
    class cellMesh;
    class faceMesh;
    class pointMesh;

    class runtimeMesh : public geometryMesh {
    private:
        /*!\brief The region list.*/
        mutable sharedPtrList<markRegion> regions_;

        using regionMapType = std::map<std::string, markRegion *>;
        mutable regionMapType regionMap_;

        void createRegions(const controller &cont) const;

        void clearRegion() noexcept;

    public:
        hur_nodiscard static const OpenHurricane::runtimeMesh &nullObject() {
            return NullRefObj::nullRef<runtimeMesh>();
        }

        inline runtimeMesh(const object &ob) : geometryMesh(ob), regions_(), regionMap_() {
            checkAndRportMeshQuality();
        }

        inline runtimeMesh(object &&ob) noexcept
            : geometryMesh(std::move(ob)), regions_(), regionMap_() {
            checkAndRportMeshQuality();
        }

        /*!\brief Contruct from decomposing mesh.*/
        runtimeMesh(const object &ob, const controller &cac);

        /*!\brief Contruct from decomposing mesh.*/
        runtimeMesh(object &&ob, const controller &cac) noexcept;

        /*!\brief Contruct from decomposing mesh.*/
        runtimeMesh(object &&ob, const controller &cac, const std::string &meshStr) noexcept;

        /*!\brief Destructor.*/
        inline virtual ~runtimeMesh() noexcept { clearRegion(); }

        /*!\breif Return const access to the cell volume.*/
        hur_nodiscard inline const realArray &cellVol() const;

        /*!\breif Return const access to the cell centre.*/
        hur_nodiscard inline const vectorArray &cellCntr() const;

        /*!\breif Return const access to the face area.*/
        hur_nodiscard inline const vectorArray &fA() const;

        /*!\breif Return const access to the face centre.*/
        hur_nodiscard inline const vectorArray &faceCntr() const;

        /*!\breif Return const access to the face weight.*/
        hur_nodiscard inline const realArray &faceWgt() const;

        /*!\brief Return const access to the region list.*/
        hur_nodiscard inline const sharedPtrList<markRegion> &regions() const noexcept;

        hur_nodiscard inline const markRegion &region(const std::string &regionName) const;
    };

} // namespace OpenHurricane

#include "runtimeMesh.inl"