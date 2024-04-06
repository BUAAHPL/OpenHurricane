#include "formFaceFromCell.hpp"
/*!
 * \file formFaceFromCell.inl
 * \brief In-Line subroutines of the <i>formFaceFromCell.hpp</i> file.
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

inline OpenHurricane::tmpFace::tmpFace() : faceIndex_(-1), leftCell_(-1), rightCell_(-1) {}

inline OpenHurricane::tmpFace::tmpFace(integer faceIndex, integer leftCell, integer rightCell)
    : faceIndex_(faceIndex), leftCell_(leftCell), rightCell_(rightCell) {}

inline OpenHurricane::tmpFace::~tmpFace() noexcept {}

inline OpenHurricane::formFaceFromCell::formFaceFromCell(cellList &cells, faceList &faces,
                                                     const cellZoneList &cellZones,
                                                     faceZoneList &faceZones,
                                                     const integerListListList &faceZoneEleConn,
                                                     const integer faceTableCapacity,
                                                     const integer tmpIntSetSize)
    : faceTable_(faceTableCapacity), tmpIntSet_(tmpIntSetSize, -1), faceId_(), cells_(cells),
      faces_(faces), cellZones_(cellZones), faceZones_(faceZones),
      faceZoneEleConn_(faceZoneEleConn) {}

inline OpenHurricane::formFaceFromCell::~formFaceFromCell() noexcept {
    clear();
}

#endif // USES_CGNS