/*!
 * \file formFaceFromCell.hpp
 * \brief Headers of forming faces from given cells.
 *        The subroutines and functions are in the <i>formFaceFromCell.cpp</i> file.
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
#ifdef USES_CGNS
#include "commonInclude.hpp"
#include "geometryMesh.hpp"
#include "integerArray.hpp"
#include "uiListHashMap.hpp"

namespace OpenHurricane {

    class tmpFace {
    public:
        integer faceIndex_;
        integer leftCell_;
        integer rightCell_;

    public:

        /**
         * \brief Null constructor.
         */
        inline tmpFace();

        inline tmpFace(integer faceIndex, integer leftCell, integer rightCell);

        inline ~tmpFace() noexcept;
    };

    class formFaceFromCell {
    private:
        intSetHashTable<tmpFace *> faceTable_;
        integerList tmpIntSet_;

        integerList faceId_;

        cellList &cells_;
        faceList &faces_;
        const cellZoneList &cellZones_;
        faceZoneList &faceZones_;

        const integerListListList &faceZoneEleConn_;

        /**
         * \brief Return the corner nodes forming the face of the tetrahedral cells.
         * \param[in] celli - The global index of cell
         * \param[in] facei - The local index of face from cell[celli]
         */
        static void tetraFace(const cell &celli, const integer facei, integerList &faceNode);

        /**
         * \brief Return the corner nodes forming the face of the pyramid cells.
         * \param[in] celli - The global index of cell
         * \param[in] facei - The local index of face from cell[celli]
         */
        static void pyraFace(const cell &celli, const integer facei, integerList &faceNode);

        /**
         * \brief Return the corner nodes forming the face of the pentahedral cells.
         * \param[in] celli - The global index of cell
         * \param[in] facei - The local index of face from cell[celli]
         */
        static void pentaFace(const cell &celli, const integer facei, integerList &faceNode);

        /**
         * \brief Return the corner nodes forming the face of the hehxahedral cells.
         * \param[in] celli - The global index of cell
         * \param[in] facei - The local index of face from cell[celli]
         */
        static void hexaFace(const cell &celli, const integer facei, integerList &faceNode);

        void formingFace();

        void reorderIntFace(const cellZone &cz, faceZone &fz, integer &offset);

        void reorderBndFace(const integerListList &eleConn, faceZone &fz, integer &offset);

        void settingFaceNodes();

    public:

        formFaceFromCell() = delete;

        inline formFaceFromCell(cellList &cells, faceList &faces, const cellZoneList &cellZones,
                                faceZoneList &faceZones, const integerListListList &faceZoneEleConn,
                                const integer faceTableCapacity, const integer tmpIntSetSize = 8);

        /**
         * \brief Destructor.
         */
        inline ~formFaceFromCell() noexcept;

        void clear() noexcept;

        void operator()();
    };
} // namespace OpenHurricane

#include "formFaceFromCell.inl"
#endif // USES_CGNS