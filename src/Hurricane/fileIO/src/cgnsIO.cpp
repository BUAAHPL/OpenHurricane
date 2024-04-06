/*!
 * \file cgnsIO.cpp
 * \brief Main subroutines of the <i>cgnsIO.hpp</i> file.
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
#include "cgnsIO.hpp"

#ifdef USES_CGNS

void OpenHurricane::cgnsIO::setMode(const int flag) {
    switch (flag) {
    case READ:
        mode_ = READ;
        break;
    case WRITE:
        mode_ = WRITE;
        break;
    case MODIFY:
        mode_ = MODIFY;
        break;
    default:
        mode_ = UNKNOWN_MODE;
        break;
    }
}

OpenHurricane::cgnsIO::CGNSGridLocation
OpenHurricane::cgnsIO::convertGridLocation(const GridLocation_t &gl) const {
    switch (gl) {
    case GridLocationNull:
        return CGNSGridLocation::CGNS_NULL;
    case GridLocationUserDefined:
        return CGNSGridLocation::CGNS_UserDefined;
    case Vertex:
        return CGNSGridLocation::CGNS_Vertex;
    case CellCenter:
        return CGNSGridLocation::CGNS_CellCenter;
    case FaceCenter:
        return CGNSGridLocation::CGNS_FaceCenter;
    case IFaceCenter:
        return CGNSGridLocation::CGNS_IFaceCenter;
    case JFaceCenter:
        return CGNSGridLocation::CGNS_JFaceCenter;
    case KFaceCenter:
        return CGNSGridLocation::CGNS_KFaceCenter;
    case EdgeCenter:
        return CGNSGridLocation::CGNS_EdgeCenter;
    default:
        LFatal("Unknown grid location type");
    }
    return CGNSGridLocation::CGNS_NULL;
}

GridLocation_t OpenHurricane::cgnsIO::convertGridLocation(const CGNSGridLocation &gl) const {
    switch (gl) {
    case CGNSGridLocation::CGNS_NULL:
        return GridLocationNull;
    case CGNSGridLocation::CGNS_UserDefined:
        return GridLocationUserDefined;
    case CGNSGridLocation::CGNS_Vertex:
        return Vertex;
    case CGNSGridLocation::CGNS_CellCenter:
        return CellCenter;
    case CGNSGridLocation::CGNS_FaceCenter:
        return FaceCenter;
    case CGNSGridLocation::CGNS_IFaceCenter:
        return IFaceCenter;
    case CGNSGridLocation::CGNS_JFaceCenter:
        return JFaceCenter;
    case CGNSGridLocation::CGNS_KFaceCenter:
        return KFaceCenter;
    case CGNSGridLocation::CGNS_EdgeCenter:
        return EdgeCenter;
    default:
        LFatal("Unknown grid location type");
    }
    return GridLocationNull;
}

OpenHurricane::cgnsIO::cgnsIO(const fileName &fN, const int flg)
    : filename_(fN), cgnsFN_(-1), mode_(MODIFY), openClosed_(CLOSED) {
    setMode(flg);
}

void OpenHurricane::cgnsIO::open() {
    if (filename_.empty()) {
        LFatal("Attempt to open file in CGNS format with null filename");
    }

    if (mode_ == UNKNOWN_MODE) {
        LFatal("Attempt to open file in CGNS format using unknown mode");
    }

    if (opened()) {
        LFatal("Attempt to open file: %s which was already opened", filename_.c_str());
    }

    const auto result = cg_open(filename_.c_str(), mode_, &cgnsFN_);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    openClosed_ = OPENED;
}

hur_nodiscard OpenHurricane::cellShapeType::shapeTypes
OpenHurricane::cgnsIO::checkElementType(const ElementType_t &type) const noexcept {
    switch (type) {
    case TETRA_4:
        return cellShapeType::shapeTypes::tetrahedral;
    case PYRA_5:
        return cellShapeType::shapeTypes::pyramid;
    case PENTA_6:
        return cellShapeType::shapeTypes::wedge;
    case HEXA_8:
        return cellShapeType::shapeTypes::hexahedral;
    case TRI_3:
        return cellShapeType::shapeTypes::triangular;
    case QUAD_4:
        return cellShapeType::shapeTypes::quadrilateral;
    case MIXED:
        return cellShapeType::shapeTypes::mixed;
    case NGON_n:
        return cellShapeType::shapeTypes ::polyhedral;
    default:
        LFatal("Unsupported element type: %s", ElementTypeName[type]);
    }
    return cellShapeType::shapeTypes::invalid;
}

hur_nodiscard ElementType_t
OpenHurricane::cgnsIO::convertElementType(const cellShapeType::shapeTypes &type) const noexcept {
    switch (type) {
    case cellShapeType::shapeTypes::mixed:
        return ElementType_t::MIXED;
        break;
    case cellShapeType::shapeTypes::triangular:
        return ElementType_t::TRI_3;
        break;
    case cellShapeType::shapeTypes::tetrahedral:
        return ElementType_t::TETRA_4;
        break;
    case cellShapeType::shapeTypes::quadrilateral:
        return ElementType_t::QUAD_4;
        break;
    case cellShapeType::shapeTypes::hexahedral:
        return ElementType_t::HEXA_8;
        break;
    case cellShapeType::shapeTypes::pyramid:
        return ElementType_t::PYRA_5;
        break;
    case cellShapeType::shapeTypes::wedge:
        return ElementType_t::PENTA_6;
        break;
        break;
    case cellShapeType::shapeTypes::polyhedral:
        return ElementType_t::NGON_n;
        break;
    case cellShapeType::shapeTypes::invalid:
    default:
        break;
    }
    return ElementType_t::ElementTypeNull;
}

hur_nodiscard OpenHurricane::faceShapeType::shapeTypes
OpenHurricane::cgnsIO::checkFaceElementType(const ElementType_t &type) const noexcept {
    switch (type) {
    case BAR_2:
        return faceShapeType::shapeTypes::linear;
    case TRI_3:
        return faceShapeType::shapeTypes::triangular;
    case QUAD_4:
        return faceShapeType::shapeTypes::quadrilateral;
    case NFACE_n:
        return faceShapeType::shapeTypes::polygonal;
    case MIXED:
        return faceShapeType::shapeTypes::mixed;
    default:
        LFatal("Unsupported face element type: %s", ElementTypeName[type]);
    }
    return faceShapeType::shapeTypes::invalid;
}

hur_nodiscard ElementType_t
OpenHurricane::cgnsIO::convertFaceElementType(const faceShapeType::shapeTypes &type) const noexcept {
    switch (type) {
    case faceShapeType::shapeTypes::mixed:
        return ElementType_t::MIXED;
        break;
    case faceShapeType::shapeTypes::linear:
        return ElementType_t::BAR_2;
        break;
    case faceShapeType::shapeTypes::triangular:
        return ElementType_t::TRI_3;
        break;
    case faceShapeType::shapeTypes::quadrilateral:
        return ElementType_t::QUAD_4;
        break;
    case faceShapeType::shapeTypes::polygonal:
        return ElementType_t::NFACE_n;
        break;
    case faceShapeType::shapeTypes::invalid:
    default:
        break;
    }
    return ElementType_t::ElementTypeNull;
}

hur_nodiscard OpenHurricane::faceBCType::bcTypes
OpenHurricane::cgnsIO::checkBCType(const BCType_t &type) const noexcept {
    switch (type) {
    case BCTypeNull:
        return faceBCType::bcTypes::WALL;
    case BCTypeUserDefined:
        return faceBCType::bcTypes::WALL;
    case BCAxisymmetricWedge:
        return faceBCType::bcTypes::AXIS;
    case BCDegenerateLine:
        return faceBCType::bcTypes::WALL;
    case BCDegeneratePoint:
        return faceBCType::bcTypes::WALL;
    case BCDirichlet:
        return faceBCType::bcTypes::PRESSUREFARFIELD;
    case BCExtrapolate:
        return faceBCType::bcTypes::WALL;
    case BCFarfield:
        return faceBCType::bcTypes::PRESSUREFARFIELD;
    case BCGeneral:
        return faceBCType::bcTypes::INTERIOR;
    case BCInflow:
        return faceBCType::bcTypes::VELOCITYINLET;
    case BCInflowSubsonic:
        return faceBCType::bcTypes::PRESSUREINLET;
    case BCInflowSupersonic:
        return faceBCType::bcTypes::PRESSUREFARFIELD;
    case BCNeumann:
        return faceBCType::bcTypes::PRESSUREFARFIELD;
    case BCOutflow:
        return faceBCType::bcTypes::VELOCITYINLET;
    case BCOutflowSubsonic:
        return faceBCType::bcTypes::PRESSUREOUTLET;
    case BCOutflowSupersonic:
        return faceBCType::bcTypes::OUTFLOW;
    case BCSymmetryPlane:
        return faceBCType::bcTypes::SYMMETRY;
    case BCSymmetryPolar:
        return faceBCType::bcTypes::AXIS;
    case BCTunnelInflow:
        return faceBCType::bcTypes::VELOCITYINLET;
    case BCTunnelOutflow:
        return faceBCType::bcTypes::OUTFLOW;
    case BCWall:
    case BCWallInviscid:
    case BCWallViscous:
    case BCWallViscousHeatFlux:
    case BCWallViscousIsothermal:
        return faceBCType::bcTypes::WALL;
    case FamilySpecified:
        return faceBCType::bcTypes::WALL;
    default:
        break;
    }
    return faceBCType::bcTypes::WALL;
}

hur_nodiscard BCType_t
OpenHurricane::cgnsIO::convertBCType(const faceBCType::bcTypes &type) const noexcept {
    switch (type) {
    case faceBCType::bcTypes::INTERIOR:
        return BCType_t::BCGeneral;
        break;
    case faceBCType::bcTypes::WALL:
        return BCType_t::BCWall;
        break;
    case faceBCType::bcTypes::PRESSUREINLET:
        return BCType_t::BCInflowSubsonic;
        break;
    case faceBCType::bcTypes::INLETVENT:
        return BCType_t::BCInflow;
        break;
    case faceBCType::bcTypes::INLETFAN:
        return BCType_t::BCInflow;
        break;
    case faceBCType::bcTypes::PRESSUREOUTLET:
        return BCType_t::BCOutflowSubsonic;
        break;
    case faceBCType::bcTypes::EXHAUSTFAN:
        return BCType_t::BCOutflowSubsonic;
        break;
    case faceBCType::bcTypes::OUTLETVENT:
        return BCType_t::BCOutflowSubsonic;
        break;
    case faceBCType::bcTypes::SYMMETRY:
        return BCType_t::BCSymmetryPlane;
        break;
    case faceBCType::bcTypes::PERIODICSHADOW:
        return BCType_t::BCSymmetryPlane;
        break;
    case faceBCType::bcTypes::PRESSUREFARFIELD:
        return BCType_t::BCFarfield;
        break;
    case faceBCType::bcTypes::VELOCITYINLET:
        return BCType_t::BCInflow;
        break;
    case faceBCType::bcTypes::PERIODIC:
        return BCType_t::BCSymmetryPlane;
        break;
    case faceBCType::bcTypes::FAN:
        return BCType_t::BCInflow;
        break;
    case faceBCType::bcTypes::POROUSJUMP:
        return BCType_t::BCGeneral;
        break;
    case faceBCType::bcTypes::RADIATOR:
        return BCType_t::BCInflow;
        break;
    case faceBCType::bcTypes::MASSFLOWINLET:
        return BCType_t ::BCInflow;
        break;
    case faceBCType::bcTypes::DETONATIONINLET:
        return BCType_t::BCInflowSupersonic;
        break;
    case faceBCType::bcTypes::OUTFLOW:
        return BCType_t::BCOutflowSupersonic;
        break;
    case faceBCType::bcTypes::AXIS:
        return BCType_t::BCSymmetryPolar;
        break;
    default:
        return BCType_t::BCGeneral;
        break;
    }
}

hur_nodiscard int OpenHurricane::cgnsIO::numNodes(const ElementType_t &ty) const {
    int npe;
    const auto result = cg_npe(ty, &npe);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return npe;
}

hur_nodiscard bool OpenHurricane::cgnsIO::isUnstructuredZone(const int B, const int Z) const {
    ZoneType_t zt;
    readZoneType(B, Z, &zt);
    return zt == Unstructured;
}

hur_nodiscard OpenHurricane::integerListList
OpenHurricane::cgnsIO::parsingElementConn(const cgsize_t *hur_restrict ele, const integer eleSize,
                                      const integer eleDataSize, const ElementType_t &tpptr) const {
    integerListList eleConn(eleSize);
    auto f = [&eleConn, ele](const int n, const integer iEle, const integer iEleDat) {
        eleConn[iEle].resize(n);
        for (integer i = 0; i < n; ++i) {
            eleConn[iEle][i] = (integer)ele[iEleDat + i] - 1;
        }
        return iEleDat + n;
    };
    integer iEleDat;
    integer npe;
    switch (tpptr) {
    case TETRA_4:
    case PYRA_5:
    case PENTA_6:
    case HEXA_8:
    case TRI_3:
    case QUAD_4:
        npe = numNodes(tpptr);
        iEleDat = 0;
        for (integer iEle = 0; iEle < eleConn.size(); ++iEle) {
            iEleDat = f(npe, iEle, iEleDat);
        }
        break;
    case MIXED:
        iEleDat = 0;
        for (integer iEle = 0; iEle < eleConn.size(); ++iEle) {
            npe = numNodes(ElementType_t(ele[iEleDat++]));
            iEleDat = f(npe, iEle, iEleDat);
        }
        break;
    default:
        errorAbortStr(("Unsupported element type: " + std::string(ElementTypeName[tpptr])));
    }
    return eleConn;
}

void OpenHurricane::cgnsIO::parsingCellElementConn(const cgsize_t *hur_restrict ele,
                                               const integer start, const integer end,
                                               const integer eleDataSize,
                                               const ElementType_t &tpptr, cellList &cells) const {
    auto f = [&cells, ele](const int n, const integer icell, const integer iEleDat) {
        cells[icell].nodesList().resize(n);
        for (integer i = 0; i < n; ++i) {
            cells[icell].nodesList()[i] = (integer)ele[iEleDat + i] - 1;
        }
        return iEleDat + n;
    };
    integer iEleDat;
    integer npe;
    switch (tpptr) {
    case TETRA_4:
    case PYRA_5:
    case PENTA_6:
    case HEXA_8:
    case TRI_3:
    case QUAD_4:
        npe = numNodes(tpptr);
        iEleDat = 0;
        for (integer icell = start; icell <= end; ++icell) {
            iEleDat = f(npe, icell, iEleDat);
            cells[icell].setType(checkElementType(tpptr));
        }
        break;
    case MIXED:
        iEleDat = 0;
        for (integer icell = start; icell <= end; ++icell) {
            cells[icell].setType(checkElementType(ElementType_t(ele[iEleDat])));
            npe = numNodes(ElementType_t(ele[iEleDat++]));
            iEleDat = f(npe, icell, iEleDat);
        }
        break;
    default:
        errorAbortStr(("Unsupported element type: " + std::string(ElementTypeName[tpptr])));
    }
}

#endif // USES_CGNS