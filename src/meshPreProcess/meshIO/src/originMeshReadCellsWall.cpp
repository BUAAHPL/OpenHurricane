/*!
 * \file originMeshReadCellsWall.cpp
 * \brief Main subroutines for computing cell wall flag.
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

#include "logFile.hpp"
#include "originMeshRead.hpp"
// compute near wall flag for cell:
// 1) iflag = 0, interior cell
// 2) iflag = 1, wall boundary cell with one node on the wall
// 3) iflag = 2, wall boundary cell with one face on the wall
void OpenHurricane::originMeshRead::computingCellWall(const integer globalNodes_,
                                                      const integer globalCells_,
                                                      const faceZoneList &faceZones_,
                                                      const faceList &faces_, cellList &cells_) {
    Pout << "    Info: computing near wall flag for each cell... " << std::endl;

    List<cellWallFlag::wallFlags> color(globalNodes_, cellWallFlag::wallFlags::interior);

    for (integer fzi = 0; fzi < faceZones_.size(); fzi++) {
        if (faceZones_[fzi].isWall()) {
            for (integer walli = faceZones_[fzi].firstIndex(); walli <= faceZones_[fzi].lastIndex();
                 walli++) {
                for (integer nodei = 0; nodei < faces_[walli].size(); nodei++) {
                    integer k = faces_[walli][nodei];
                    color[k] = cellWallFlag::wallFlags::oneNodeOnWall;
                }
            }
        }
    }

    for (integer celli = 0; celli < globalCells_; celli++) {
        for (integer nodei = 0; nodei < cells_[celli].nodeSize(); nodei++) {
            integer k = cells_[celli].nodei(nodei);
            if (color[k] == cellWallFlag::wallFlags::oneNodeOnWall) {
                cells_[celli].setWallFlag(color[k]);
                break;
            }
        }
    }

    for (integer fzi = 0; fzi < faceZones_.size(); fzi++) {
        if (faceZones_[fzi].isWall()) {
            for (integer walli = faceZones_[fzi].firstIndex(); walli <= faceZones_[fzi].lastIndex();
                 walli++) {
                integer celli = faces_[walli].leftCell();
                cells_[celli].setWallFlag(cellWallFlag::wallFlags::oneFaceOnWall);
            }
        }
    }
}
