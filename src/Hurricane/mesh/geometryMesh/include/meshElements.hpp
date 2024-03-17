/*!
 * \file meshElements.hpp
 * \brief Header of mesh elements, i.e., point, cell and face.
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
#include "Lists.hpp"
#include "geometryModel.hpp"
#include "vectorArray.hpp"

namespace OpenHurricane {
    /**
     * \brief The point's definition is a three dimensional vector.
     */
    using point = vector;
    /**
     * \brief The point field.
     */
    using pointField = vectorArray;

    class cell;

    /**
     * \brief The list of the cell.
     */
    using cellList = List<cell>;

    class face;

    using faceList = List<face>;

    namespace cellWallFlag {
        enum wallFlags : short {
            interior,      // Interior cell
            oneNodeOnWall, // Wall boundary cell with one node on the wall
            oneFaceOnWall  // Wall boundary cell with one face on the wall
        };

        hur_nodiscard wallFlags getWallFlags(const short wallFlag) noexcept;

    } // namespace cellWallFlag

    namespace cellShapeType {
        enum shapeTypes : short {
            mixed = 0,
            triangular = 1,
            tetrahedral = 2,
            quadrilateral = 3,
            hexahedral = 4,
            pyramid = 5,
            wedge = 6,
            polyhedral = 7,
            invalid
        };

        hur_nodiscard shapeTypes getshapeType(const short shapeType) noexcept;
    } // namespace cellShapeType

    namespace cellFacePair {

        static const std::map<short, short> HexCellFacePairMap =
            createMap<short, short>(short(0), short(1))(short(1), short(0))(short(2), short(3))(
                short(3), short(2))(short(4), short(5))(short(5), short(4));

        static const std::map<short, short> QuaCellFacePairMap = createMap<short, short>(
            short(0), short(2))(short(2), short(0))(short(1), short(3))(short(3), short(1));

    } //  namespace cellFacePair

    /**
     * \brief The class of the cell.
     */
    class cell {
    private:
        /**\brief The node lists forming the cell.*/
        integerList nodes_;

        /**\brief The face lists forming the cell.*/
        integerList faces_;

        /**\brief Original index when reading file.*/
        integer originIndex_;

        /**\brief The wall flag of cell.*/
        cellWallFlag::wallFlags wallFlag_;

        /**\brief The geometry bcType the cell.*/
        cellShapeType::shapeTypes shapeType_;

    public:
        inline cell()
            : nodes_(), faces_(), originIndex_(), wallFlag_(cellWallFlag::wallFlags::interior),
              shapeType_(cellShapeType::shapeTypes::hexahedral) {}

        cell(const short shapeType);

        cell(const cell &);

        inline ~cell() noexcept {}

        /**\brief Return the cell points.*/
        hur_nodiscard OpenHurricane::pointField points(const pointField meshPoint) const;

        hur_nodiscard inline const integerList &nodesList() const noexcept { return nodes_; }
        hur_nodiscard inline integerList &nodesList() noexcept { return nodes_; }

        /**\brief Return the integer of the ith node that this cell contains .*/
        hur_nodiscard inline integer nodei(const integer i) const noexcept {
            return nodes_.operator[](i);
        }
        hur_nodiscard inline integer nodeSize() const noexcept { return nodes_.size(); }

        hur_nodiscard inline const integerList &facesList() const noexcept { return faces_; }
        hur_nodiscard inline integerList &facesList() noexcept { return faces_; }

        /**\brief Return the integer of the ith face that this cell contains.*/
        hur_nodiscard hur_nodiscard inline integer facei(const integer i) const noexcept {
            return faces_.operator[](i);
        }

        hur_nodiscard inline integer faceSize() const noexcept { return faces_.size(); }

        hur_nodiscard inline cellWallFlag::wallFlags wallFlag() const noexcept { return wallFlag_; }

        hur_nodiscard inline cellShapeType::shapeTypes shapeType() const noexcept {
            return shapeType_;
        }

        hur_nodiscard inline integer orignIndex() const noexcept { return originIndex_; }

        hur_nodiscard integer findFace(const integer faceI) const;

        /*!\brief Return the opposite face index*/
        hur_nodiscard integer oppositeFace(const integer);

        /**\brief set nodes size and face size.*/
        inline void resize(const short NN, const short NF) {
            nodes_.resize(NN);
            faces_.resize(NF);
        }

        /**\brief set nodes size and face size.*/
        void resize();

        cell &operator=(const cell &);
        void transfer(cell &c) noexcept;

        inline void setType(const cellShapeType::shapeTypes ct) noexcept { shapeType_ = ct; }

        inline void setType(const short shapeType) noexcept {
            shapeType_ = cellShapeType::getshapeType(shapeType);
        }

        inline void setWallFlag(const cellWallFlag::wallFlags wallFlag) noexcept {
            wallFlag_ = wallFlag;
        }

        inline void setWallFlag(const short wallFlag) noexcept {
            wallFlag_ = cellWallFlag::getWallFlags(wallFlag);
        }

        inline void setOrignIndex(const integer id) noexcept { originIndex_ = id; }

        /**
         * \brief Return the corresponding canonical positions of the face in this cell.
         * \param[in] faacei - The face index in the mesh.
         */
        hur_nodiscard integer faceConoId(const integer facei) const noexcept;
        
        /**
         * \brief Return the corresponding canonical positions of the face in this cell under the CGNS format.
         * \param[in] faacei - The face index in the mesh.
         */
        hur_nodiscard integer cgnsFaceConoId(const integer facei) const noexcept;

        hur_nodiscard inline bool isTriangular() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::triangular);
        }

        hur_nodiscard inline bool isTetrahedral() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::tetrahedral);
        }

        hur_nodiscard inline bool isQuadrilateral() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::quadrilateral);
        }

        hur_nodiscard inline bool isHexahedral() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::hexahedral);
        }

        hur_nodiscard inline bool isPyramid() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::pyramid);
        }

        hur_nodiscard inline bool isWedge() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::wedge);
        }

        hur_nodiscard inline bool isPolyhedral() const noexcept {
            return (shapeType_ == cellShapeType::shapeTypes::polyhedral);
        }
    };

    // Face boundary type.
    namespace faceBCType {

        /*!\brief The bcType of boundary condition.*/
        enum bcTypes : short {
            CUTFACE = 1,          // parallel cut face
            INTERIOR = 2,         // interior
            WALL = 3,             // wall
            PRESSUREINLET = 4,    // pressure inlet
            INLETVENT = 41,       // inlet vent
            INLETFAN = 42,        // intake fan
            PRESSUREOUTLET = 5,   // pressure outlet
            EXHAUSTFAN = 51,      // exhaust fan
            OUTLETVENT = 52,      // outlet vent
            SYMMETRY = 7,         // symmetry
            PERIODICSHADOW = 8,   // periodic shadow
            PRESSUREFARFIELD = 9, // pressure far field
            VELOCITYINLET = 10,   // velocity inlet
            PERIODIC = 12,        // periodic
            FAN = 14,             // fan
            POROUSJUMP = 141,     // porous jump
            RADIATOR = 142,       // radiator
            MASSFLOWINLET = 20,   // mass flow inlet
            DETONATIONINLET = 21, //detonation inlet
            INTERFACE = 24,       // inetrface
            PARENT = 31,          // parent (hanging node)
            OUTFLOW = 36,         // outflow
            AXIS = 37             // axis
        };

        hur_nodiscard bcTypes getBcTypes(const short bcType) noexcept;

        static std::map<std::string, bcTypes> bcTypeStringMap =
            createMap<std::string, bcTypes>("CUTFACE", CUTFACE)("INTERIOR", INTERIOR)("WALL", WALL)(
                "PRESSUREINLET", PRESSUREINLET)("INLETVENT", INLETVENT)("INLETFAN", INLETFAN)(
                "PRESSUREOUTLET", PRESSUREOUTLET)("EXHAUSTFAN", EXHAUSTFAN)(
                "OUTLETVENT", OUTLETVENT)("SYMMETRY", SYMMETRY)("PERIODICSHADOW", PERIODICSHADOW)(
                "PRESSUREFARFIELD", PRESSUREFARFIELD)("VELOCITYINLET", VELOCITYINLET)(
                "PERIODIC", PERIODIC)("FAN", FAN)("POROUSJUMP", POROUSJUMP)("RADIATOR", RADIATOR)(
                "MASSFLOWINLET", MASSFLOWINLET)("DETONATIONINLET", DETONATIONINLET)(
                "INTERFACE", INTERFACE)("PARENT", PARENT)("OUTFLOW", OUTFLOW)("AXIS", AXIS);
    } // namespace faceBCType

    // Face shape type.
    namespace faceShapeType {
        enum shapeTypes : short {
            mixed = 0,
            linear = 2,
            triangular = 3,
            quadrilateral = 4,
            polygonal = 5,
            invalid
        };

        hur_nodiscard shapeTypes getShapeTypes(const short shapeType) noexcept;
    } // namespace faceShapeType

    /**
     * \brief The class of the face.
     */
    class face : public integerList {
    public:
        using Base = integerList;

    private:
        /**\brief The ID of the left cell.*/
        integer leftCell_;

        /**\brief The ID of the right cell.*/
        integer rightCell_;

        /**\brief The boundary bcType of the face.*/
        faceBCType::bcTypes type_;

    public:
        static void areaAndCentre(const pointField &p, vector &fA, vector &fC);

        /**\brief Construct from null.*/
        inline face()
            : integerList(), leftCell_(0), rightCell_(0), type_(faceBCType::bcTypes::INTERIOR) {}

        /**\brief Construct from components.*/
        inline face(const integer lCell, const integer rCell, const integerList &node)
            : integerList(node), leftCell_(lCell), rightCell_(rCell),
              type_(faceBCType::bcTypes::INTERIOR) {}

        /**\brief Construct from components.*/
        inline face(const integer lCell, const integer rCell, const integerList &node,
                    const faceBCType::bcTypes type)
            : integerList(node), leftCell_(lCell), rightCell_(rCell), type_(type) {}

        /**\brief Construct from components.*/
        face(const integer lCell, const integer rCell, const faceBCType::bcTypes type,
             const integer nn);

        /**\brief Construct as copy.*/
        inline face(const face &f)
            : integerList(f), leftCell_(f.leftCell_), rightCell_(f.rightCell_), type_(f.type_) {}

        /**\brief Destructor.*/
        inline ~face() noexcept {}

        /**\brief Return the const access to the ID of the left cell.*/
        hur_nodiscard inline const integer &leftCell() const noexcept { return leftCell_; }

        /**\brief Return the non-const access to the ID of the left cell.*/
        hur_nodiscard inline integer &leftCell() noexcept { return leftCell_; }

        /**\brief Return the const access to the ID of the right cell.*/
        hur_nodiscard inline const integer &rightCell() const noexcept { return rightCell_; }

        /**\brief Return the non-const access to the ID of the right cell.*/
        hur_nodiscard inline integer &rightCell() noexcept { return rightCell_; }

        /**\brief Use with care*/
        hur_forceinline hur_nodiscard integer oppositeCell(const integer cellI) const noexcept {
            return leftCell_ + rightCell_ - cellI;
        }

        /**\brief Return access to the boundary bcType of the face.*/
        hur_nodiscard inline faceBCType::bcTypes type() const noexcept { return type_; }

        /**\brief Set the boundary bcType for the face.*/
        inline void setBCType(const short bT) noexcept { type_ = faceBCType::getBcTypes(bT); }

        /**\brief Return true if this face zone is parellel cut face zone.*/
        hur_nodiscard inline bool isCutFace() const noexcept {
            return (type_ == faceBCType::bcTypes::CUTFACE);
        }

        /**\brief Return true if this face is interior face
         * The interface face is considered as interior face.*/
        hur_nodiscard inline bool isInterior() const noexcept {
            return ((type_ == faceBCType::bcTypes::INTERIOR) ||
                    (type_ == faceBCType::bcTypes::INTERFACE));
        }

        /**\brief Return true if this face zone is interface face zone.*/
        hur_nodiscard inline bool isInterface() const noexcept {
            return (type_ == faceBCType::bcTypes::INTERFACE);
        }

        /**\brief Return true if this face zone is wall face zone.*/
        hur_nodiscard inline bool isWall() const noexcept {
            return (type_ == faceBCType::bcTypes::WALL);
        }

        /**\brief Return true if this face zone is periodic face zone.*/
        hur_nodiscard inline bool isPeriodic() const noexcept {
            return (type_ == faceBCType::bcTypes::PERIODIC);
        }

        /**\brief Return true if this face zone is periodic-shadow face zone.*/
        hur_nodiscard inline bool isPeriodicShadow() const noexcept {
            return (type_ == faceBCType::bcTypes::PERIODICSHADOW);
        }

        hur_nodiscard inline bool isLinear() const noexcept { return (this->size() == 2); }
        hur_nodiscard inline bool isTriangular() const noexcept { return (this->size() == 3); }
        hur_nodiscard inline bool isQuadrilateral() const noexcept { return (this->size() == 4); }
        hur_nodiscard inline bool isPolygonal() const noexcept { return (this->size() > 4); }

        /*!\brief Return if this face is connected with the given one by an edge*/
        hur_nodiscard bool coEdge(const face &) const;

        face &operator=(const face &f);

        void transfer(face &a) noexcept;
    };

    class faceCalc {
    public:
        static inline point triangularArea(const point &point1, const point &point2,
                                           const point &point3) {
            return (real(0.5) * ((point2 - point1) ^ (point3 - point1)));
        }

        static inline point triangularCentre(const point &point1, const point &point2,
                                             const point &point3) {
            return (real(1.0 / 3.0) * (point1 + point2 + point3));
        }

        static inline point quadrilateralArea(const point &point1, const point &point2,
                                              const point &point3, const point &point4) {
            return (real(0.5) * ((point3 - point1) ^ (point2 - point4)));
        }

        static inline point quadrilateralCentre(const point &point1, const point &point2,
                                                const point &point3, const point &point4) {
            return (real(0.25) * (point1 + point2 + point3 + point4));
        }

        static inline point projectionInface(const point &spacePoint, const point &faceCentre,
                                             const vector &faceNormal) {
            const auto nn = faceNormal.normalized();
            const auto u = spacePoint - faceCentre;
            const auto ppu = u - (u * nn) * nn;
            return ppu + faceCentre;
        }

        static inline bool pointInTriangularFace(const point &pointP, const point &point1,
                                                 const point &point2, const point &point3) {
            const auto PA = point1 - pointP;
            const auto PB = point2 - pointP;
            const auto PC = point3 - pointP;

            const auto t1 = PA ^ PB;
            const auto t2 = PB ^ PC;
            const auto t3 = PC ^ PA;
            return (t1 * t2 >= 0) && (t1 * t3 >= 0) && (t2 * t3 >= 0);
        }

        static inline bool pointInQuadrilateralFace(const point &pointP, const point &point1,
                                                    const point &point2, const point &point3,
                                                    const point &point4) {
            const auto PA = point1 - pointP;
            const auto PB = point2 - pointP;
            const auto PC = point3 - pointP;
            const auto PD = point4 - pointP;

            const auto t1 = PA ^ PB;
            const auto t2 = PB ^ PC;
            const auto t3 = PC ^ PD;
            const auto t4 = PD ^ PA;
            return (t1 * t2 >= 0) && (t1 * t3 >= 0) && (t1 * t4 >= 0) && (t2 * t3 >= 0) &&
                   (t2 * t4 >= 0) && (t3 * t4 >= 0);
        }
    };

    /**
     * \brief Return the minimum integer of the nodes that the face contains.
     * \param[in] f The face.
     * \return The minimum integer of the nodes
     */
    hur_nodiscard inline integer min(const face &f) {
        integer i = f[0];
        for (int faceI = 1; faceI < f.size(); faceI++) {
            i = min(f[faceI], i);
        }
        return i;
    }

    /**
     * \brief Return the maximum integer of the nodes that the face contains.
     * \param[in] f The face.
     * \return The maximum integer of the nodes
     */
    hur_nodiscard inline integer max(const face &f) {
        integer i = f[0];
        for (int faceI = 1; faceI < f.size(); faceI++) {
            i = max(f[faceI], i);
        }
        return i;
    }

    /*!\brief Template specialization for feature<face>.*/
    template <> class feature<face> {
        face p_;

    public:
        declareClassNames;

        /*!\brief Construct from primitive.*/
        explicit feature(const face &p) : p_(p) {}
    };
} // namespace OpenHurricane
