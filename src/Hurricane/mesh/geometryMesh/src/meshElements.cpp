/*!
 * \file meshElements.cpp
 * \brief The subroutines and functions of  mesh elements, i.e., point, cell and face.
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
#include "meshElements.hpp"

hur_nodiscard OpenHurricane::cellWallFlag::wallFlags
OpenHurricane::cellWallFlag::getWallFlags(const short wallFlag) noexcept {
    switch (wallFlag) {
    case cellWallFlag::interior:
        return cellWallFlag::wallFlags::interior;
        break;
    case cellWallFlag::oneNodeOnWall:
        return cellWallFlag::wallFlags::oneNodeOnWall;
        break;
    case cellWallFlag::oneFaceOnWall:
        return cellWallFlag::wallFlags::oneFaceOnWall;
        break;
    default:
        return cellWallFlag::wallFlags::interior;
        break;
    }
}

hur_nodiscard OpenHurricane::cellShapeType::shapeTypes
OpenHurricane::cellShapeType::getshapeType(const short shapeType) noexcept {
    switch (shapeType) {
    case cellShapeType::mixed:
        return cellShapeType::shapeTypes::mixed;
        break;
    case cellShapeType::triangular:
        return cellShapeType::shapeTypes::triangular;
        break;
    case cellShapeType::tetrahedral:
        return cellShapeType::shapeTypes::tetrahedral;
        break;
    case cellShapeType::quadrilateral:
        return cellShapeType::shapeTypes::quadrilateral;
        break;
    case cellShapeType::hexahedral:
        return cellShapeType::shapeTypes::hexahedral;
        break;
    case cellShapeType::pyramid:
        return cellShapeType::shapeTypes::pyramid;
        break;
    case cellShapeType::wedge:
        return cellShapeType::shapeTypes::wedge;
        break;
    case cellShapeType::polyhedral:
        return cellShapeType::shapeTypes::polyhedral;
        break;
    default:
        return cellShapeType::shapeTypes::invalid;
        break;
    }
}

OpenHurricane::cell::cell(const short shapeType)
    : nodes_(), faces_(), originIndex_(), wallFlag_(cellWallFlag::wallFlags::interior) {
    setType(shapeType);
}

OpenHurricane::cell::cell(const cell &c)
    : nodes_(c.nodes_), faces_(c.faces_), originIndex_(c.originIndex_), wallFlag_(c.wallFlag_),
      shapeType_(c.shapeType_) {}

hur_nodiscard inline OpenHurricane::pointField
OpenHurricane::cell::points(const pointField meshPoints) const {
    pointField p(nodes_.size());

    for (integer ni = 0; ni < nodes_.size(); ni++) {
        p[ni] = meshPoints[nodei(ni)];
    }
    return p;
}

//set nodes size and face size
void OpenHurricane::cell::resize() {
    if (isTriangular()) {
        resize(3, 3);
    } else if (isTetrahedral()) {
        resize(4, 4);
    } else if (isQuadrilateral()) {
        resize(4, 4);
    } else if (isHexahedral()) {
        resize(8, 6);
    } else if (isPyramid()) {
        resize(5, 5);
    } else if (isWedge()) {
        resize(6, 5);
    } else if (isPolyhedral()) {
        LFatal("Must specify the number of the nodes and faces for plolyhedral cell");
    }
}

hur_nodiscard OpenHurricane::integer OpenHurricane::cell::findFace(const integer faceI) const {
    integer cellFaceIndex = -1;
    for (integer localI = 0; localI < this->facesList().size(); ++localI) {
        if (faceI == this->facesList()[localI]) {
            cellFaceIndex = localI;
        }
    }

    return cellFaceIndex;
}

hur_nodiscard OpenHurricane::integer OpenHurricane::cell::oppositeFace(const integer fI) {
    if ((this->shapeType_ != cellShapeType::shapeTypes::hexahedral) &&
        (this->shapeType_ != cellShapeType::shapeTypes::quadrilateral)) {
        LFatal("Not support cell type %d", this->shapeType_);
        return -1;
    } else {
        integer cellFaceIndex = findFace(fI);

        std::map<short, short>::const_iterator iter;
        iter = cellFacePair::HexCellFacePairMap.find(short(cellFaceIndex));

        if (iter == cellFacePair::HexCellFacePairMap.end()) {
            LFatal("Not find opposite face, maybe cell type error ");
        }

        return this->facesList()[integer(iter->second)];
    }
}
OpenHurricane::cell &OpenHurricane::cell::operator=(const cell &c) {
    if (this != std::addressof(c)) {
        nodes_.clear();
        nodes_.resize(c.nodes_.size());
        nodes_ = c.nodes_;

        faces_.clear();
        faces_.resize(c.faces_.size());
        faces_ = c.faces_;

        wallFlag_ = c.wallFlag_;
        shapeType_ = c.shapeType_;
        originIndex_ = c.originIndex_;
    }
    return *this;
}

void OpenHurricane::cell::transfer(cell &c) noexcept {
    if (this == std::addressof(c)) {
        return;
    }
    nodes_.transfer(c.nodes_);
    faces_.transfer(c.faces_);
    wallFlag_ = c.wallFlag_;
    shapeType_ = c.shapeType_;
    originIndex_ = c.originIndex_;
}

hur_nodiscard OpenHurricane::integer OpenHurricane::cell::faceConoId(const integer facei) const noexcept {
    for (integer i = 0; i < faces_.size(); ++i) {
        if (faces_[i] == facei) {
            return i;
        }
    }
    return -1;
}

hur_nodiscard OpenHurricane::integer
OpenHurricane::cell::cgnsFaceConoId(const integer facei) const noexcept {
    const auto ii = faceConoId(facei);
    if (ii == -1) {
        return -1;
    }
    if (isHexahedral()) {
        if (ii == 0) {
            return 4;
        } else if (ii == 1) {
            return 2;
        } else if (ii == 2) {
            return 1;
        } else if (ii == 3) {
            return 3;
        } else if (ii == 4) {
            return 0;
        } else if (ii == 5) {
            return 5;
        }
    } else if (isTetrahedral()) {
        if (ii == 0) {
            return 0;
        } else if (ii == 1) {
            return 1;
        } else if (ii == 2) {
            return 3;
        } else if (ii == 3) {
            return 2;
        }
    } else if (isPyramid()) {
        if (ii == 0) {
            return 4;
        } else if (ii == 1) {
            return 2;
        } else if (ii == 2) {
            return 1;
        } else if (ii == 3) {
            return 3;
        } else if (ii == 4) {
            return 0;
        }
    } else if (isWedge()) {
        if (ii == 0) {
            return 2;
        } else if (ii == 1) {
            return 1;
        } else if (ii == 2) {
            return 0;
        } else if (ii == 3) {
            return 3;
        } else if (ii == 4) {
            return 4;
        }
    }
    return ii;
}

hur_nodiscard OpenHurricane::faceBCType::bcTypes
OpenHurricane::faceBCType::getBcTypes(const short bcType) noexcept {
    switch (bcType) {
    case bcTypes::CUTFACE:
        return bcTypes::CUTFACE;
        break;
    case bcTypes::INTERIOR:
        return bcTypes::INTERIOR;
        break;
    case bcTypes::WALL:
        return bcTypes::WALL;
        break;
    case bcTypes::PRESSUREINLET:
        return bcTypes::PRESSUREINLET;
        break;
    case bcTypes::INLETVENT:
        return bcTypes::INLETVENT;
        break;
    case bcTypes::INLETFAN:
        return bcTypes::INLETFAN;
        break;
    case bcTypes::PRESSUREOUTLET:
        return bcTypes::PRESSUREOUTLET;
        break;
    case bcTypes::EXHAUSTFAN:
        return bcTypes::EXHAUSTFAN;
        break;
    case bcTypes::OUTLETVENT:
        return bcTypes::OUTLETVENT;
        break;
    case bcTypes::SYMMETRY:
        return bcTypes::SYMMETRY;
        break;
    case bcTypes::PERIODICSHADOW:
        return bcTypes::PERIODICSHADOW;
        break;
    case bcTypes::PRESSUREFARFIELD:
        return bcTypes::PRESSUREFARFIELD;
        break;
    case bcTypes::VELOCITYINLET:
        return bcTypes::VELOCITYINLET;
        break;
    case bcTypes::PERIODIC:
        return bcTypes::PERIODIC;
        break;
    case bcTypes::FAN:
        return bcTypes::FAN;
        break;
    case bcTypes::POROUSJUMP:
        return bcTypes::POROUSJUMP;
        break;
    case bcTypes::RADIATOR:
        return bcTypes::RADIATOR;
        break;
    case bcTypes::MASSFLOWINLET:
        return bcTypes::MASSFLOWINLET;
        break;
    case bcTypes::DETONATIONINLET:
        return bcTypes::DETONATIONINLET;
        break;
    case bcTypes::INTERFACE:
        return bcTypes::INTERFACE;
        break;
    case bcTypes::PARENT:
        return bcTypes::PARENT;
        break;
    case bcTypes::OUTFLOW:
        return bcTypes::OUTFLOW;
        break;
    case bcTypes::AXIS:
        return bcTypes::AXIS;
        break;
    default:
        return bcTypes::WALL;
        break;
    }
}
hur_nodiscard OpenHurricane::faceShapeType::shapeTypes
OpenHurricane::faceShapeType::getShapeTypes(const short shapeType) noexcept {
    switch (shapeType) {
    case shapeTypes::mixed:
        return shapeTypes::mixed;
        break;
    case shapeTypes::linear:
        return shapeTypes::linear;
        break;
    case shapeTypes::triangular:
        return shapeTypes::triangular;
        break;
    case shapeTypes::quadrilateral:
        return shapeTypes::quadrilateral;
        break;
    case shapeTypes::polygonal:
        return shapeTypes::polygonal;
        break;
    case shapeTypes::invalid:
        return shapeTypes::invalid;
        break;
    default:
        return shapeTypes::invalid;
        break;
    }
}

void OpenHurricane::face::areaAndCentre(const pointField &p, vector &fA, vector &fC) {
    integer nFPoints = p.size();

    if (nFPoints == 3) {
        fC = (real(1.0 / 3.0)) * (p[0] + p[1] + p[2]);
        fA = (real(0.5)) * ((p[1] - p[0]) ^ (p[2] - p[0]));
    } else if (nFPoints == 4) {
        fC = (real(0.25)) * (p[0] + p[1] + p[2] + p[3]);
        fA = (real(0.5) * ((p[2] - p[0]) ^ (p[3] - p[1])));
    } else {
        vector sumN = Zero;
        real sumA = real(0.0);
        vector sumAc = Zero;

        point fCentre = p[0];
        for (integer pi = 1; pi < nFPoints; pi++) {
            fCentre += p[pi];
        }
        fCentre /= real(nFPoints);

        for (integer pi = 0; pi < nFPoints; pi++) {
            const point &nextPoint = p[(pi + 1) % nFPoints];

            vector c = p[pi] + nextPoint + fCentre;
            vector n = (nextPoint - p[pi]) ^ (fCentre - p[pi]);
            real a = mag(n);

            sumN += n;
            sumA += a;
            sumAc += a * c;
        }

        if (sumA < rootVeryTiny) {
            fC = fCentre;
            fA = Zero;
        } else {
            fC = real(1.0 / 3.0) * sumAc / sumA;
            fA = real(0.5) * sumN;
        }
    }
}

OpenHurricane::face::face(const integer lCell, const integer rCell, const faceBCType::bcTypes Type,
                      const integer nn)
    : leftCell_(lCell), rightCell_(rCell), type_(Type) {

    resize(nn);
}

hur_nodiscard bool OpenHurricane::face::coEdge(const face &f) const {
    integer publicPoint = 0;
    for (integer i = 0; i < this->size(); ++i) {
        for (integer j = 0; j < f.size(); ++j) {
            if (this->operator[](i) == f[j]) {
                publicPoint++;
            }
        }
    }
    return publicPoint == 2;
}

OpenHurricane::face &OpenHurricane::face::operator=(const face &f) {
    if (this != std::addressof(f)) {
        Base::resize(f.size());
        this->leftCell_ = f.leftCell_;
        this->rightCell_ = f.rightCell_;
        Base::operator=(f);
        this->type_ = f.type_;
    }
    return *this;
}

void OpenHurricane::face::transfer(face &f) noexcept {
    if (this != std::addressof(f)) {
        return;
    }
    integerList::transfer(f);
    this->leftCell_ = f.leftCell_;
    this->rightCell_ = f.rightCell_;
    this->type_ = f.type_;
}
