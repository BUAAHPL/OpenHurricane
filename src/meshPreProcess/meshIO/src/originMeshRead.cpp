/*!
 * \file originMeshRead.cpp
 * \brief Main subroutines for reading origin mesh.
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

#include "originMeshRead.hpp"
#include <sstream>

namespace OpenHurricane {
    createClassNameStr(originMeshRead, "originMeshRead");
    createObjFty(originMeshRead, controller);
} // namespace OpenHurricane

OpenHurricane::originMeshRead::originMeshRead()
    : meshStr_(), fileName_(), nProcs_(1), ND_(0), origin_(Zero), axis_(Zero), globalNodes_(0),
      globalFaces_(0), globalInterior_(0), globalCells_(0), points_(), faces_(), cells_(),
      pointZones_(), faceZones_(), cellZones_(), periodicPairZone_(), hasBeenRead_(false),
      hasPeriodic_(false), hasHangingNodes_(false), hasInterface_(false), originCellIndex_(),
      interiorWallFaceMap_(), decomposeList_(), originMeshDecompSize_(1), cellLoadWeights_() {}

OpenHurricane::originMeshRead::originMeshRead(const fileName &fN, const int nP)
    : meshStr_(), fileName_(fN), nProcs_(nP), ND_(0), origin_(Zero), axis_(Zero), globalNodes_(0),
      globalFaces_(0), globalInterior_(0), globalCells_(0), points_(), faces_(), cells_(),
      pointZones_(), faceZones_(), cellZones_(), periodicPairZone_(), hasBeenRead_(false),
      hasPeriodic_(false), hasHangingNodes_(false), hasInterface_(false), originCellIndex_(),
      interiorWallFaceMap_(), decomposeList_(), originMeshDecompSize_(1), cellLoadWeights_() {}

OpenHurricane::originMeshRead::originMeshRead(const std::string &str, const int nP)
    : meshStr_(str), fileName_(), nProcs_(nP), ND_(0), origin_(Zero), axis_(Zero), globalNodes_(0),
      globalFaces_(0), globalInterior_(0), globalCells_(0), points_(), faces_(), cells_(),
      pointZones_(), faceZones_(), cellZones_(), periodicPairZone_(), hasBeenRead_(false),
      hasPeriodic_(false), hasHangingNodes_(false), hasInterface_(false), originCellIndex_(),
      interiorWallFaceMap_(), decomposeList_(), originMeshDecompSize_(1), cellLoadWeights_() {}

OpenHurricane::uniquePtr<OpenHurricane::originMeshRead>
OpenHurricane::originMeshRead::creator(const fileName &fN, const int nP, const controller &cont) {
    string meshFormat = cont.findWord("meshFormat");

    defineInObjCreator(originMeshRead, static_cast<std::string>(meshFormat), controller, (fN, nP));
}

void OpenHurricane::originMeshRead::clear() noexcept {
    meshStr_.clear();
    fileName_.clear();
    points_.clear();
    faces_.clear();
    cells_.clear();

    pointZones_.clear();

    faceZones_.clear();
    cellZones_.clear();
    periodicPairZone_.clear();
    originCellIndex_.clear();
    interiorWallFaceMap_.clear();
}

void OpenHurricane::originMeshRead::read(string &gridUnit) {
    if (!fileName_.empty()) {
        Pout << "    Buffering mesh from file: " << fileName_ << "..." << std::endl;
    }
    reading(gridUnit);
}

void OpenHurricane::originMeshRead::printMeshInfo() const {
    Pout << std::endl
         << "    Mesh info:" << std::endl
         << "      Total nodes: " << globalNodes_ << std::endl
         << "      Total faces: " << globalFaces_ << std::endl
         << "      Total cells: " << globalCells_ << std::endl
         << "      Zones info:" << std::endl;
    std::streamsize mw = Pout.width();
    Pout.setf(std::ios::right);
    for (integer pzi = 0; pzi < pointZones_.size(); ++pzi) {
        Pout.width(20);
        Pout << pointZones_[pzi].name() << ": ";
        Pout.width(10);
        Pout << pointZones_[pzi].size() << " (nodes)" << std::endl;
    }
    for (integer fzi = 0; fzi < faceZones_.size(); ++fzi) {
        Pout.width(20);
        Pout << faceZones_[fzi].name() << ": ";
        Pout.width(10);
        Pout << faceZones_[fzi].size() << " (faces)" << std::endl;
    }

    for (integer czi = 0; czi < cellZones_.size(); ++czi) {
        Pout.width(20);
        Pout << cellZones_[czi].name() << ": ";
        Pout.width(10);
        Pout << cellZones_[czi].size() << " (cells)" << std::endl;
    }
    Pout.unsetf(std::ios::right);
    Pout << std::endl;
}