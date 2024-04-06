/*!
 * \file formFaceFromCell.cpp
 * \brief Main subroutines for forming faces from given cells.
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

#include "formFaceFromCell.hpp"
#ifdef USES_CGNS

void OpenHurricane::formFaceFromCell::tetraFace(const cell &celli, const integer facei,
                                                integerList &faceNode) {
    faceNode.resize(3);
    auto f = [&celli, &faceNode](const integer n1, const integer n2, const integer n3) {
        faceNode[0] = celli.nodei(n1);
        faceNode[1] = celli.nodei(n2);
        faceNode[2] = celli.nodei(n3);
    };

    switch (facei) {
    case 0: // Face F1: N0, N2, N1
        f(0, 2, 1);
        break;
    case 1: // Face F2: N0, N1, N3
        f(0, 1, 3);
        break;
    case 2: // Face F3: N1, N2, N3
        f(1, 2, 3);
        break;
    case 3: // Face F4: N2, N0, N3
        f(2, 0, 3);
        break;
    default:
        errorAbortStr(("Invalid local face index: " + toString(facei) +
                       ". The local face index for a tetrahedral cell can only be 0~3"));
        break;
    }
}

void OpenHurricane::formFaceFromCell::pyraFace(const cell &celli, const integer facei,
                                               integerList &faceNode) {
    auto f = [&celli, &faceNode](const integer n1, const integer n2, const integer n3) {
        faceNode.resize(3);
        faceNode[0] = celli.nodei(n1);
        faceNode[1] = celli.nodei(n2);
        faceNode[2] = celli.nodei(n3);
    };

    switch (facei) {
    case 0: // Face F1: N0, N3, N2, N1
        faceNode.resize(4);
        faceNode[0] = celli.nodei(0);
        faceNode[1] = celli.nodei(3);
        faceNode[2] = celli.nodei(2);
        faceNode[3] = celli.nodei(1);
        break;
    case 1: // Face F2: N0, N1, N4
        f(0, 1, 4);
        break;
    case 2: // Face F3: N1, N2, N4
        f(1, 2, 4);
        break;
    case 3: // Face F4: N2, N3, N4
        f(2, 3, 4);
        break;
    case 4: // Face F5: N3, N0, N4
        f(3, 0, 4);
        break;
    default:
        errorAbortStr(("Invalid local face index: " + toString(facei) +
                       ". The local face index for a pyramid cell can only be 0~4"));
        break;
    }
}

void OpenHurricane::formFaceFromCell::pentaFace(const cell &celli, const integer facei,
                                                integerList &faceNode) {
    auto f3 = [&celli, &faceNode](const integer n1, const integer n2, const integer n3) {
        faceNode.resize(3);
        faceNode[0] = celli.nodei(n1);
        faceNode[1] = celli.nodei(n2);
        faceNode[2] = celli.nodei(n3);
    };

    auto f4 = [&celli, &faceNode](const integer n1, const integer n2, const integer n3,
                                  const integer n4) {
        faceNode.resize(4);
        faceNode[0] = celli.nodei(n1);
        faceNode[1] = celli.nodei(n2);
        faceNode[2] = celli.nodei(n3);
        faceNode[3] = celli.nodei(n4);
    };

    switch (facei) {
    case 0: // Face F1: N0, N1, N4, N3
        f4(0, 1, 4, 3);
        break;
    case 1: // Face F2: N1, N2, N5, N4
        f4(1, 2, 5, 4);
        break;
    case 2: // Face F3: N2, N0, N3, N5
        f4(2, 0, 3, 5);
        break;
    case 3: // Face F4: N0, N2, N1
        f3(0, 2, 1);
        break;
    case 4: // Face F5: N3, N4, N5
        f3(3, 4, 5);
        break;
    default:
        errorAbortStr(("Invalid local face index: " + toString(facei) +
                       ". The local face index for a pentahedral cell can only be 0~4"));
        break;
    }
}

void OpenHurricane::formFaceFromCell::hexaFace(const cell &celli, const integer facei,
                                               integerList &faceNode) {
    faceNode.resize(4);
    auto f = [&celli, &faceNode](const integer n1, const integer n2, const integer n3,
                                 const integer n4) {
        faceNode[0] = celli.nodei(n1);
        faceNode[1] = celli.nodei(n2);
        faceNode[2] = celli.nodei(n3);
        faceNode[3] = celli.nodei(n4);
    };

    switch (facei) {
    case 0: // Face F1: N0, N3, N2, N1
        f(0, 3, 2, 1);
        break;
    case 1: // Face F2: N0, N1, N5, N4
        f(0, 1, 5, 4);
        break;
    case 2: // Face F3: N1, N2, N6, N5
        f(1, 2, 6, 5);
        break;
    case 3: // Face F4: N2, N3, N7, N6
        f(2, 3, 7, 6);
        break;
    case 4: // Face F5: N0, N4, N7, N3
        f(0, 4, 7, 3);
        break;
    case 5: // Face F6: N4, N5, N6, N7
        f(4, 5, 6, 7);
        break;
    default:
        errorAbortStr(("Invalid local face index: " + toString(facei) +
                       ". The local face index for a hexahedral cell can only be 0~5"));
        break;
    }
}

void OpenHurricane::formFaceFromCell::formingFace() {
    integer fcount = 0;
    auto func = [this, &fcount](const integer fsize, const integer ci,
                                void (*faceSetFunc)(const cell &, const integer, integerList &)) {
        cells_[ci].facesList().resize(fsize);
        for (integer fi = 0; fi < fsize; ++fi) {
            (*faceSetFunc)(cells_[ci], fi, tmpIntSet_);
            auto iter = faceTable_.find(tmpIntSet_);
            if (iter == faceTable_.end()) {
                cells_[ci].facesList()[fi] = fcount;
                faceTable_.emplace(tmpIntSet_, new tmpFace(fcount, ci, -1));
                fcount++;
            } else {
                cells_[ci].facesList()[fi] = iter->second->faceIndex_;
                iter->second->rightCell_ = ci;
            }
        }
    };
    integer fzi = 0;
    integer offset = 0;
    for (const auto &e : cellZones_) {
        for (integer i = e.firstIndex(); i <= e.lastIndex(); ++i) {
            switch (cells_[i].shapeType()) {
            case cellShapeType::shapeTypes::tetrahedral:
                func(4, i, tetraFace);
                break;
            case cellShapeType::shapeTypes::pyramid:
                func(5, i, pyraFace);
                break;
            case cellShapeType::shapeTypes::wedge:
                func(5, i, pentaFace);
                break;
            case cellShapeType::shapeTypes::hexahedral:
                func(6, i, hexaFace);
                break;
            case cellShapeType::shapeTypes::triangular:
            case cellShapeType::shapeTypes::quadrilateral:
            default:
                LFatal("Unsupported cell type");
                break;
            }
        }
        if (faceZones_[fzi].name() != "int_" + e.name()) {
            errorAbortStr(("The face zone: " + faceZones_[fzi].name() +
                           " is not interior face zone of cell zone: " + e.name()));
        }

        fzi++;
    }

    faces_.resize(fcount);
    faceId_.resize(fcount, -1);
}

void OpenHurricane::formFaceFromCell::reorderIntFace(const cellZone &cz, faceZone &fz,
                                                     integer &offset) {
    integer interiorCount = 0;
    for (auto &e : faceTable_) {
        if (e.second->leftCell_ != -1 && e.second->rightCell_ != -1) {
            if (e.second->leftCell_ >= cz.firstIndex() && e.second->leftCell_ <= cz.lastIndex() &&
                e.second->rightCell_ >= cz.firstIndex() && e.second->rightCell_ <= cz.lastIndex()) {
                const auto fid = offset + interiorCount;
                faceId_[e.second->faceIndex_] = fid;
                interiorCount++;

                faces_[fid].leftCell() = e.second->leftCell_;
                faces_[fid].rightCell() = e.second->rightCell_;
            }
        }
    }
    fz.setFirstIndex(offset);
    fz.setLastIndex(offset + interiorCount - 1);
    offset += interiorCount;
}

void OpenHurricane::formFaceFromCell::reorderBndFace(const integerListList &eleConn, faceZone &fz,
                                                     integer &offset) {
    fz.setFirstIndex(offset);
    fz.setLastIndex(offset + eleConn.size() - 1);
    integer fcount = 0;
    for (const auto &ec : eleConn) {
        auto iter = faceTable_.find(ec);
        if (iter != faceTable_.end()) {
            const auto fid = offset + fcount;
            faceId_[iter->second->faceIndex_] = fid;
            faces_[fid].leftCell() = iter->second->leftCell_;
            faces_[fid].rightCell() = iter->second->rightCell_;

            ++fcount;
        } else {
            LFatal("Face missed in table");
        }
    }
    offset += eleConn.size();
}

void OpenHurricane::formFaceFromCell::settingFaceNodes() {
    auto func1 = [](const cell &celli, const integer fn, integerList &ndl,
                    void (*faceSetFunc)(const cell &, const integer, integerList &)) {
        (*faceSetFunc)(celli, fn, ndl);
    };

    auto func2 = [this, func1](const cell &celli, const integer fn) {
        switch (celli.shapeType()) {
        case cellShapeType::shapeTypes::tetrahedral:
            func1(celli, fn, tmpIntSet_, tetraFace);
            break;
        case cellShapeType::shapeTypes::pyramid:
            func1(celli, fn, tmpIntSet_, pyraFace);
            break;
        case cellShapeType::shapeTypes::wedge:
            func1(celli, fn, tmpIntSet_, pentaFace);
            break;
        case cellShapeType::shapeTypes::hexahedral:
            func1(celli, fn, tmpIntSet_, hexaFace);
            break;
        case cellShapeType::shapeTypes::triangular:
        case cellShapeType::shapeTypes::quadrilateral:
        default:
            LFatal("Unsupported cell type");
            break;
        }
    };
    tmpIntSet_ = -1;
    integer ci = 0;
    boolList faceSet(faceId_.size(), false);
    for (auto &e : cells_) {
        integer fn = 0;
        for (auto &fi : e.facesList()) {
            const auto fid = faceId_[fi];
            if (!faceSet[fi]) {
                func2(e, fn);
                faces_[fid].resize(tmpIntSet_.size());

                if (ci == faces_[fid].leftCell()) {
                    for (integer j = 0; j < faces_[fid].size(); ++j) {
                        faces_[fid][j] = tmpIntSet_[faces_[fid].size() - j - 1];
                    }
                } else if (ci == faces_[fid].rightCell()) {
                    for (integer j = 0; j < faces_[fid].size(); ++j) {
                        faces_[fid][j] = tmpIntSet_[j];
                    }
                } else {
                    LFatal("Face not right");
                }
                faceSet[fi] = true;
            }
            fi = fid;
            fn++;
        }
        ci++;
    }
}

void OpenHurricane::formFaceFromCell::clear() noexcept {
    for (auto &e : faceTable_) {
        HurDelete(e.second);
    }
    faceTable_.clear();
}

void OpenHurricane::formFaceFromCell::operator()() {
    formingFace();
    integer fzi = 0;
    integer offset = 0;
    for (const auto &e : cellZones_) {
        reorderIntFace(e, faceZones_[fzi], offset);
        fzi++;
    }

    for (integer i = fzi; i < faceZones_.size(); ++i) {
        reorderBndFace(faceZoneEleConn_[i - fzi], faceZones_[i], offset);
    }
    settingFaceNodes();
}
#endif // USES_CGNS