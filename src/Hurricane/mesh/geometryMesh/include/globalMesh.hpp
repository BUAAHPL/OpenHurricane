/*!
 * \file globalMesh.hpp
 * \brief Headers of the global mesh.
 *        The subroutines and functions are in the <i>globalMesh.cpp</i> file.
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
#include "decomposingByMetis.hpp"
#include "geometryMesh.hpp"
#include "globalIndeces.hpp"
#include "smartPointer.hpp"

namespace OpenHurricane {

    class globalMesh {
    private:
        const geometryMesh &mesh_;

        mutable globalCellIndex *globalCellIndecesPtr_;

        mutable globalFaceIndex *globalFaceIndecesPtr_;

        mutable globalPointIndex *globalPointIndecesPtr_;

    public:

        /*!\brief Construct from components.*/
        inline globalMesh(const geometryMesh &gMesh);

        inline globalMesh(const geometryMesh &gMesh, globalCellIndex &gCI, globalFaceIndex &gFI,
                          globalPointIndex &gPI);

        /**\brief Construct from components.*/
        inline globalMesh(const geometryMesh &gMesh, const decomposingMeshByMetis &dcM);

        // Destructor
        ~globalMesh() noexcept;

        /*!\brief Return const access to the mseh.*/
        hur_nodiscard inline const geometryMesh &mesh() const noexcept;

        /*!\brief Clear.*/
        void clear() noexcept;

        /**\brief Return global cell index.*/
        hur_nodiscard inline integer toGlobalCellIndex(const integer zoneid,
                                                       const integer localIndex) const;

        /**\brief Return global cell index.*/
        hur_nodiscard inline integer toGlobalCellIndex(const integer localIndex) const;

        /**\brief Return global face index.*/
        hur_nodiscard inline integer toGlobalFaceIndex(const integer zoneId,
                                                       const integer localId) const;

        /**\brief Return global point index.*/
        hur_nodiscard inline integer toGlobalPointIndex(const integer zoneId,
                                                        const integer localId) const;

        /**\brief Return global point index.*/
        hur_nodiscard inline integer toGlobalPointIndex(const integer localId) const;

        hur_nodiscard inline const globalCellIndex &globalCellIndeces() const noexcept;

        hur_nodiscard inline const globalFaceIndex &globalFaceIndeces() const noexcept;

        hur_nodiscard inline const globalPointIndex &globalPointIndeces() const noexcept;

        void gatherFaceZone(const integer zonei, faceZone &gfz, integerArrayArray &faceConnect,
                            integer &_offset) const;
    };

} // namespace OpenHurricane

#include "globalMesh.inl"