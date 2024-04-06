/*!
 * \file originMeshRead.hpp
 * \brief Headers of the origin mesh reading.
 *        The subroutines and functions are in the <i>originMeshRead.cpp</i> file.
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

#include "commonInclude.hpp"
#include "geometryMesh.hpp"
#include "integerArray.hpp"
#include "iteration.hpp"

namespace OpenHurricane {

    class hdf5I;

    class originMeshRead {
    protected:
        /**\brief Mesh file string.*/
        std::string meshStr_;

        /*!\brief The filename of the mesh file.*/
        fileName fileName_;

        /*!\brief The numberof all processors.*/
        int nProcs_;

        /*!\brief The dimension set of the origin mesh.*/
        unsigned short int ND_;

        vector origin_;
        vector axis_;

        /*!\brief The number of the mesh nodes.*/
        integer globalNodes_;

        /*!\brief The number of the mesh faces.*/
        integer globalFaces_;

        /*!\brief The number of the interior faces.*/
        integer globalInterior_;

        /*!\brief The number of the mesh cells.*/
        integer globalCells_;

        /*!\brief The global mesh points.*/
        pointField points_;

        /*!\brief The global mesh faces.*/
        faceList faces_;

        /*!\brief The global mesh cells.*/
        cellList cells_;

        /*!\brief Zone information.*/

        /*!\brief The global point zones.*/
        pointZoneList pointZones_;

        /*!\brief The global face zones.*/
        faceZoneList faceZones_;

        /*!\brief The global cell zones.*/
        cellZoneList cellZones_;

        /*!\brief The global periodic zone.*/
        periodicPairList periodicPairZone_;

        /*!\brief The mesh file has been read.*/
        bool hasBeenRead_;

        /*!\brief If the grid has the periodic boundary?.*/
        bool hasPeriodic_;

        /*!\brief If the grid has the hanging nodes?.*/
        bool hasHangingNodes_;

        bool hasInterface_;

        /*!\brief Origin cell index of .h5 file mesh.*/
        integerArrayArray originCellIndex_;

        std::map<integer, integer> interiorWallFaceMap_;

        /**\brief The lists of the processors used in the origin mesh (the input mesh).*/
        integerArrayArray decomposeList_;
        /**\brief The number of the processors used in the origin mesh (the input mesh).*/
        integer originMeshDecompSize_;

    protected:
        /*!\brief To compute forming faces for cell.*/
        //void formingFaces();

        /*!\brief To set the node and face size of each cell.*/
        static void assignCellType(cellList &);

        /*!\brief To compute forming nodes for cell.*/
        //void formingNodes();
        static void triangularCellNodes(const integer celli, const faceList &, cellList &);
        static void tetrahedralCellNodes(const integer, const faceList &, cellList &);
        static void quadrilateralCellNodes(const integer celli, const faceList &, cellList &);
        static void pyramidCellNodes(const integer, const faceList &, cellList &);
        static void wedgeCellNodes(const integer, const faceList &, cellList &);
        static void hexahedralCellNodes(const integer, const faceList &, cellList &);
        static void polyhedralCellNodes(const integer celli, const faceList &, cellList &);

        static bool isContained(const integerList &l1, const integerList &l2);

    public:
        declareClassNames;
        declareObjFty(originMeshRead, controller, (const fileName &fN, const int nP), (fN, nP));

        /*!\brief Null constructor.*/
        originMeshRead();

        /*!\brief Construct from components.*/
        explicit originMeshRead(const fileName &, const int);

        originMeshRead(const std::string &str, const int);

        /*!\brief Select null constructed.*/
        static uniquePtr<originMeshRead> creator(const fileName &, const int,
                                                 const controller &cont);

        /**\brief Destructor.*/
        virtual ~originMeshRead() noexcept { clear(); }

        inline const vector &origin() const noexcept;

        inline const vector &axis() const noexcept;

        inline const fileName &meshFileName() const;

        inline int numberOfProcessor() const;

        inline short meshDimensionSet() const;

        inline integer globalNodes() const;

        inline integer globalFaces() const;

        inline integer globalInterior() const;

        inline integer globalCells() const;

        inline const pointField &points() const;

        inline const faceList &faces() const;

        inline const cellList &cells() const;

        inline const pointZoneList &pointZones() const;

        inline const faceZoneList &faceZones() const;

        inline const cellZoneList &cellZones() const;

        inline const periodicPairList &periodicPairZones() const;

        inline integerArrayArray &originCellIndex();

        inline bool hasBeenRead() const;

        /*!\brief Return true if the grid has periodic boundary.*/
        inline bool hasPeriodic() const;

        /*!\brief Return true if the grid has hanging nodes.*/
        inline bool hasHangingNodes() const;

        inline bool hasInterface() const;

        static void formingFaces(const integer, const faceList &, cellList &);

        static void formingNodes(const integer, const faceList &, cellList &);

        /*!\brief To compute cell wall flag.*/
        static void computingCellWall(const integer globalNodes_, const integer globalCells_,
                                      const faceZoneList &faceZones_, const faceList &faces_,
                                      cellList &cells_);

        virtual void clear() noexcept;

        inline std::map<integer, integer> &interiorWallFaceMap();

        inline const std::map<integer, integer> &interiorWallFaceMap() const;

    public:
        void read(string &gridUnit);

        inline virtual hur_nodiscard bool checkCoupledWall() const noexcept;
        inline virtual hur_nodiscard bool isHurricaneMesh() const noexcept;

    protected:
        virtual void reading(string &gridUnit) {}

        void printMeshInfo() const;

    public:
        /**\brief The lists of the processors used in the origin mesh (the input mesh).*/
        inline const integerArrayArray &decomposeList() const;
        /**\brief The number of the processors used in the origin mesh (the input mesh).*/
        inline integer originMeshDecompSize() const noexcept;

    protected:
        /**
         * \brief Cell load weights.
         */
        integerArray cellLoadWeights_;

    public:
        /**
         * \brief Cell load weights.
         */
        hur_nodiscard inline const integerArray &cellLoadWeights() const noexcept;
    };

} // namespace OpenHurricane

#include "originMeshRead.inl"