/*!
 * \file basemMesh.hpp
 * \brief Header of base mesh parameters.
 *       The subroutines and functions are in the <i>baseMesh.cpp</i> file.
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

#include "controller.hpp"
#include "dataStructure.hpp"
#include "integer.hpp"
#include "meshElements.hpp"
#include "realArray.hpp"
#include "vectorArray.hpp"
#include "zoneMesh.hpp"

namespace OpenHurricane {

    /**
     * \brief Mesh checking warning threshold.
     */
    class baseMeshWarningThreshold {
    public:
        /*!\brief Aspect ratio warning threshold.*/
        static constexpr real aspectRatioThreshold = real(10000);

        /*!\brief Volume ratio warning threshold.*/
        static constexpr real volumeRatioThreshold = real(30);

        /*!\brief Face skewness warning threshold.*/
        static constexpr real skewenessThreshold = real(4);

        /*!\brief Face Orthogonality warning threshold.*/
        static real OrthogonalityThreshold;
    };

    class baseMeshCoreNums {
    protected:
        /**\brief Number of points.*/
        integer nPoints_;

        /**\brief Number of ghost points.*/
        integer nGhostPoints_;

        /**\brief Number of internal points (<=nPoints_).*/
        integer nInternalPoints_;

        /**\brief Number of faces.*/
        integer nFaces_;

        /**\brief Number of interior faces (<=nFace_).*/
        integer nInteriorFaces_;

        /**\brief Number of boundary faces.*/
        integer nBoundFaces_;

        /**\brief Number of ghost faces.*/
        integer nGhostFaces_;

        /**\brief Number of cells (internal cells).*/
        integer nCells_;

        /**\brief Number of ghost cells.*/
        integer nGhostCells_;

        /**\brief Number of total cells ( = nCells_ + nGhostCells_ ).*/
        integer nTotalCells_;

    public:
        inline baseMeshCoreNums()
            : nPoints_(0), nGhostPoints_(0), nInternalPoints_(0), nFaces_(0), nInteriorFaces_(0),
              nBoundFaces_(0), nGhostFaces_(0), nCells_(0), nGhostCells_(0), nTotalCells_(0) {}

        inline baseMeshCoreNums(const integer nPoint, const integer nFace,
                                const integer nInteriorFace, const integer nCell)
            : nPoints_(nPoint), nGhostPoints_(0), nInternalPoints_(0), nFaces_(nFace),
              nInteriorFaces_(nInteriorFace), nBoundFaces_(0), nGhostFaces_(0), nCells_(nCell),
              nGhostCells_(0), nTotalCells_(0) {
            nBoundFaces_ = nFaces_ - nInteriorFaces_;
            nGhostCells_ = nBoundFaces_;
        }

        baseMeshCoreNums(const baseMeshCoreNums &) = delete;
        baseMeshCoreNums &operator=(const baseMeshCoreNums &) = delete;

        inline ~baseMeshCoreNums() noexcept {}

        inline void setGhostCellSize(const integer newSize) noexcept { nGhostCells_ = newSize; }

        inline void setTotalCellSize(const integer newSize) noexcept { nTotalCells_ = newSize; }

        // Return the number of points (except ghost points).
        hur_nodiscard inline integer nPoints() const noexcept { return nPoints_; }

        hur_nodiscard inline integer nInternalPoints() const noexcept { return nInternalPoints_; }

        // Return the number of ghost points.
        hur_nodiscard inline integer nGhostPoints() const noexcept { return nGhostPoints_; }

        // Return the total number of points.
        hur_nodiscard inline integer nTotalPoints() const noexcept {
            return (nPoints_ + nGhostPoints_);
        }

        // Return the number of faces (except ghost faces).
        hur_nodiscard inline integer nFaces() const noexcept { return nFaces_; }

        hur_nodiscard inline integer nInteriorFaces() const noexcept { return nInteriorFaces_; }

        hur_nodiscard inline integer nBoundFaces() const noexcept { return nBoundFaces_; }

        // Return the number of ghost faces.
        hur_nodiscard inline integer nGhostFaces() const noexcept { return nGhostFaces_; }

        // Return the total number of faces.
        hur_nodiscard inline integer nTotalFaces() const noexcept {
            return (nFaces_ + nGhostFaces_);
        }

        // Return the internal cells number.
        hur_nodiscard inline integer nCells() const noexcept { return nCells_; }

        // Return the ghost cells number.
        hur_nodiscard inline integer nGhostCells() const noexcept { return nGhostCells_; }

        // Return the total number of cells,
        // including internal cells and ghost cells.
        hur_nodiscard inline integer nTotalCells() const noexcept { return nTotalCells_; }

        // Return the total number of cells
        hur_nodiscard inline integer size() const noexcept { return nTotalCells_; }
    };

    class baseMeshCoreGeomInfo {
    protected:
        /**\brief Cell centres.*/
        mutable uniquePtr<vectorArray> cellCentrePtr_;
        virtual void calcCellCentreAndVol() const = 0;

        /*
         *\brief Vector pointing from the left side cell of face i to the right cell.
         * i.e. adjoinCellCtr[i] = cellCentre[cr] - cellCentre[cl]
         *		where cr is the adjoined right cell of face i
         *			  cl is the adjoined left cell of face i.
         */
        mutable uniquePtr<vectorArray> adjoinCellCtrPtr_;
        virtual void calcAdjoinnCellCentre() const = 0;

        /**\brief Face centres.*/
        mutable uniquePtr<vectorArray> faceCenterPtr_;
        virtual void calcFaceCentreAndArea() const = 0;

        /**\brief Vector from face centre to the left adjoin cell centre.
         *	i.e. faceCtrToLeftCellCtr[fi] = faceCenter[fi] - cellCentre[cl].
         *		 where cr is the adjoined right cell of face i
         *			   cl is the adjoined left cell of face i.
         */
        mutable uniquePtr<vectorArray> faceCtrToLeftCellCtrPtr_;
        virtual void calcFaceCtrToCellCtr() const = 0;

        /**\brief Vector from face centre to the right adjoin cell centre.
         *	i.e. faceCtrToRightCellCtr[fi] = faceCenter[fi] - cellCentre[cr].
         *		 where cr is the adjoined right cell of face i
         *			   cl is the adjoined left cell of face i.
         */
        mutable uniquePtr<vectorArray> faceCtrToRightCellCtrPtr_;

        /**\brief Cell volumes.*/
        mutable uniquePtr<realArray> cellVolumePtr_;

        /**\brief Faces areas.*/
        mutable uniquePtr<vectorArray> faceAreaPtr_;

        /**\brief Face weights for construct the value on the face
         *  U_ij = weight_ij * U_i + (1 - weight_ij) * U_j.*/
        mutable uniquePtr<realArray> faceWeightPtr_;

        virtual void calcFaceWeight() const = 0;

        /*!
         * \brief The node-cell weights pointer
         * (primarily for node-based green gauss gradient scheme).
         */
        mutable uniquePtr<realArrayArray> nodeWeightPtr_;

        virtual void calcNodeCellWeight() const = 0;

        /**\brief The aspect ratio of cells.*/
        mutable uniquePtr<realArray> aspectRatioPtr_;
        virtual void calcCellAspectRatio() const = 0;

        /**\brief  The grid spacing.*/
        mutable uniquePtr<realArray> deltaMaxPtr_;
        virtual void calcCellDeltaMax() const = 0;

        // Geometric shapes

        /**\brief Are all cells hexahedral in 3D or quadrilateral in 2D.*/
        bool areAllHexOrQuadri_;

        // Neighbours

        /**\brief The neighbour cells list of the cell i (only for internal cells).*/
        mutable uniquePtr<integerListList> CNCPtr_;
        virtual void calcCellNeighbourCells() const = 0;

        /**\brief The neighbour cells list of the point.*/
        mutable uniquePtr<integerListList> PNCPtr_;
        virtual void calcPointNeighbourCells() const = 0;

        short secondNeighbourCellsSize_;

        /**\brief Cell weights for gradient calcuation when leastSquareGrad applied.*/
        mutable uniquePtr<vectorArrayArray> cellWeightPtr_;
        virtual void calcCellWeight(const string &key) const = 0;

        /**\brief The origin of the coordinate system of the mesh.*/
        point origin_;

        /**\brief The end of axis of the coordinate system of the mesh.*/
        point axis_;

        /*!\brief Minimum grid size.*/
        real minSize_;

        /**\brief The second neighbour cells of the face i (only for internal faces and boundary faces, not including ghost faces)
         * (only work when all of the cells that the mesh contains are hexahedral).*/
        mutable uniquePtr<integerListList> SNCFPtr_;

        void clearMesh() noexcept {
            cellCentrePtr_.clear();
            adjoinCellCtrPtr_.clear();
            faceCenterPtr_.clear();
            faceCtrToLeftCellCtrPtr_.clear();
            faceCtrToRightCellCtrPtr_.clear();
            cellVolumePtr_.clear();
            faceAreaPtr_.clear();
            faceWeightPtr_.clear();
            nodeWeightPtr_.clear();
            aspectRatioPtr_.clear();
            CNCPtr_.clear();
            PNCPtr_.clear();
            SNCFPtr_.clear();
            cellWeightPtr_.clear();
            deltaMaxPtr_.clear();
        }

    public:
        inline baseMeshCoreGeomInfo()
            : cellCentrePtr_(nullptr), adjoinCellCtrPtr_(nullptr), faceCenterPtr_(nullptr),
              faceCtrToLeftCellCtrPtr_(nullptr), faceCtrToRightCellCtrPtr_(nullptr),
              cellVolumePtr_(nullptr), faceAreaPtr_(nullptr), faceWeightPtr_(nullptr),
              nodeWeightPtr_(nullptr), aspectRatioPtr_(nullptr), deltaMaxPtr_(nullptr),
              areAllHexOrQuadri_(false), CNCPtr_(nullptr), PNCPtr_(nullptr),
              secondNeighbourCellsSize_(short(2)), cellWeightPtr_(nullptr), origin_(Zero),
              axis_(1.0, 0.0, 0.0), minSize_(0.0), SNCFPtr_(nullptr) {}

        baseMeshCoreGeomInfo(const baseMeshCoreGeomInfo &) = delete;
        baseMeshCoreGeomInfo &operator=(const baseMeshCoreGeomInfo &) = delete;

        inline ~baseMeshCoreGeomInfo() noexcept {}

        /**\brief The origin of the coordinate system of the mesh.*/
        hur_nodiscard inline point &origin() noexcept { return origin_; }

        /**\brief The origin of the coordinate system of the mesh.*/
        hur_nodiscard inline const point &origin() const noexcept { return origin_; }

        /**\brief The end of axis of the coordinate system of the mesh.*/
        hur_nodiscard inline point &axis() noexcept { return axis_; }

        /**\brief The end of axis of the coordinate system of the mesh.*/
        hur_nodiscard inline const point &axis() const noexcept { return axis_; }

        hur_nodiscard inline const realArray &faceWeight() const {
            if (!faceWeightPtr_) {
                calcFaceWeight();
            }
            return *faceWeightPtr_;
        }

        hur_nodiscard inline const realArrayArray &NCWgt() const {
            if (!nodeWeightPtr_) {
                calcNodeCellWeight();
            }
            return *nodeWeightPtr_;
        }

        hur_nodiscard inline const realArray &aspectRatio() const {
            if (!aspectRatioPtr_) {
                calcCellAspectRatio();
            }
            return *aspectRatioPtr_;
        }

        /**\brief Return cell spacing field.*/
        hur_nodiscard inline const realArray &deltaMax() const {
            if (!deltaMaxPtr_) {
                calcCellDeltaMax();
            }
            return *deltaMaxPtr_;
        }

        /**\brief Return cell centre field.*/
        hur_nodiscard inline const vectorArray &cellCentre() const {
            if (!cellCentrePtr_) {
                calcCellCentreAndVol();
            }
            return *cellCentrePtr_;
        }

        /*
         *\brief Vector pointing from the left side cell of face i to the right cell.
         * i.e. adjoinCellCtr[i] = cellCentre[cr] - cellCentre[cl]
         *		where cr is the adjoined right cell of face i
         *			  cl is the adjoined left cell of face i.
         */
        hur_nodiscard inline const vectorArray &adjoinCellCtr() const {
            if (!adjoinCellCtrPtr_) {
                calcAdjoinnCellCentre();
            }
            return *adjoinCellCtrPtr_;
        }

        /**\brief Return face centre field.*/
        hur_nodiscard inline const vectorArray &faceCentre() const {
            if (!faceCenterPtr_) {
                calcFaceCentreAndArea();
            }
            return *faceCenterPtr_;
        }

        /**\brief Vector from face centre to the left adjoin cell centre.
         *	i.e. faceCtrToLeftCellCtr[fi] = faceCenter[fi] - cellCentre[cl].
         *		 where cr is the adjoined right cell of face i
         *			   cl is the adjoined left cell of face i.
         */
        hur_nodiscard inline const vectorArray &faceCtrToLeftCellCtr() const {
            if (!faceCtrToLeftCellCtrPtr_) {
                calcFaceCtrToCellCtr();
            }
            return *faceCtrToLeftCellCtrPtr_;
        }

        /**\brief Vector from face centre to the right adjoin cell centre.
         *	i.e. faceCtrToRightCellCtr[fi] = faceCenter[fi] - cellCentre[cr].
         *		 where cr is the adjoined right cell of face i
         *			   cl is the adjoined left cell of face i.
         */
        hur_nodiscard inline const vectorArray &faceCtrToRightCellCtr() const {
            if (!faceCtrToRightCellCtrPtr_) {
                calcFaceCtrToCellCtr();
            }
            return *faceCtrToRightCellCtrPtr_;
        }

        /**\brief Return cell volume field.*/
        hur_nodiscard inline const realArray &cellVolume() const {
            if (!cellVolumePtr_) {
                calcCellCentreAndVol();
            }
            return *cellVolumePtr_;
        }

        /**\brief Return face area field.*/
        hur_nodiscard inline const vectorArray &faceArea() const {
            if (!faceAreaPtr_) {
                calcFaceCentreAndArea();
            }
            return *faceAreaPtr_;
        }

        /**\brief Return cell weight fieldfield.*/
        hur_nodiscard inline const vectorArrayArray &cWgt(const string &key) const {
            if (!cellWeightPtr_) {
                calcCellWeight(key);
            }
            return *cellWeightPtr_;
        }

        /**\brief Return the neighbour cells list of the cell (only for internal cells).*/
        hur_nodiscard inline const integerListList &cellNeighbourCells() const {
            if (!CNCPtr_) {
                calcCellNeighbourCells();
            }
            return *CNCPtr_;
        }

        hur_nodiscard inline const integerList &cellNeighbourCells(const integer celli) const {
            if (hasCellNeighbourCells()) {
                return cellNeighbourCells()[celli];
            } else {
                calcCellNeighbourCells();
                return cellNeighbourCells()[celli];
            }
        }

        /**\brief Return the neighbour cells list of the point.*/
        hur_nodiscard inline const integerListList &pointNeighbourCells() const {
            if (!PNCPtr_) {
                calcPointNeighbourCells();
            }
            return *PNCPtr_;
        }
        hur_nodiscard const integerList &pointNeighbourCells(const integer pointi) const {
            if (hasPointNeighbourCells()) {
                return pointNeighbourCells()[pointi];
            } else {
                calcPointNeighbourCells();
                return pointNeighbourCells()[pointi];
            }
        }

        hur_nodiscard inline bool hasCellNeighbourCells() const noexcept {
            return !(CNCPtr_.isNull());
        }
        hur_nodiscard inline bool hasPointNeighbourCells() const noexcept {
            return !(PNCPtr_.isNull());
        }

        hur_nodiscard inline bool hascellCentre() const noexcept {
            return !(cellCentrePtr_.isNull());
        }

        hur_nodiscard inline bool hasfaceCentre() const noexcept {
            return !(faceCenterPtr_.isNull());
        }

        hur_nodiscard inline bool hascellVolume() const noexcept {
            return !(cellVolumePtr_.isNull());
        }

        hur_nodiscard inline bool hasfaceArea() const noexcept {
            return !(faceAreaPtr_.isNull());
        }

        hur_nodiscard inline bool hasfaceWeight() const noexcept {
            return !(faceWeightPtr_.isNull());
        }

        hur_nodiscard inline bool hasSecondNeighbourCellFaces() const noexcept {
            return !(SNCFPtr_.isNull());
        }
    };

    class baseMeshCore : public baseMeshCoreNums, public baseMeshCoreGeomInfo {
    protected:
        /*!\brief Points.*/
        pointField points_;

        /*!\brief Faces.*/
        faceList faces_;

        /*!\brief Cells.*/
        cellList cells_;

        // Zone information

        /*!\brief Point zones.*/
        pointZoneList pointZones_;

        /*!\brief Face zones.*/
        faceZoneList faceZones_;

        /*!\brief Cut face zones.*/
        cutZoneList cutZones_;

        /*!\brief periodic face zones.*/
        perZoneList perZones_;

        /*!\brief Cell zones.*/
        cellZoneList cellZones_;

    public:
        inline baseMeshCore()
            : baseMeshCoreNums(), baseMeshCoreGeomInfo(), points_(), faces_(), cells_(),
              pointZones_(), faceZones_(), cutZones_(), perZones_(), cellZones_() {}

        baseMeshCore(const baseMeshCore &) = delete;
        baseMeshCore &operator=(const baseMeshCore &) = delete;

        inline ~baseMeshCore() noexcept {}

        /**\brief Return mesh points.*/
        hur_nodiscard inline const pointField &points() const noexcept { return points_; }

        /**\brief Return faces.*/
        hur_nodiscard inline const faceList &faces() const noexcept { return faces_; }

        /**\brief Return cells.*/
        hur_nodiscard inline const cellList &cells() const noexcept { return cells_; }

        // Zone information

        /**\brief Return point zones.*/
        hur_nodiscard inline const pointZoneList &pointZones() const noexcept {
            return pointZones_;
        }

        /**\brief Return face zones.*/
        hur_nodiscard inline const faceZoneList &faceZones() const noexcept { return faceZones_; }

        /**\brief Return cut face zones.*/
        hur_nodiscard inline const cutZoneList &cutZones() const noexcept { return cutZones_; }

        /**\brief Return periodic face zones.*/
        hur_nodiscard inline const perZoneList &perZones() const noexcept { return perZones_; }

        /**\brief Return cell zones.*/
        hur_nodiscard inline const cellZoneList &cellZones() const noexcept { return cellZones_; }
    };

    /**
     * \brief The base class of base mesh.
     */
    class baseMesh : public baseMeshCore {
    public:
        enum class cellNeighbourCellType : short {
            faceNeighbour = 0,
            twoFaceNeighbour = 1,
            nodeNeighbour = 2
        };

    protected:
        /**\brief Case controller file.*/
        const controller &cont_;

    private:
        /**
         * \brief  The bcType of ghost cell.
         *        -# 0 for no ghost cells;
         *        -# 1 for one layer of ghost cell layers;
         *        -# 2 for two layers of ghost cell layers.
         *        -# 3 for three layers of ghost cell layers, and so on.
         */
        short ghostCellType_;

        bool isGhostCellCreated_;

        /*\brief Can or cannot form two or more ghost cell layers*/
        bool twoLayer_;

        cellNeighbourCellType cellNeighbourCell_;

        std::map<std::string, cellNeighbourCellType> cellNeighbourCellTypeMap_;

        void createCellNeighbourCellTypeMap();

        // Private Member Functions

        /**\brief Calculate the neighbour cells list of the cell.*/
        virtual void calcCellNeighbourCells() const;

        /**\brief Calculate the neighbour cells list of the cell.*/
        virtual void calcPointNeighbourCells() const;

    protected:

        /**\brief Calculate face centre and area.*/
        virtual void calcFaceCentreAndArea() const;

        /**\brief Vector from face centre to the adjoin cell centre.
         *	i.e. faceCtrToCellCtr[fi][0] = faceCenter[fi] - cellCentre[cl].
         *	     faceCtrToCellCtr[fi][1] = faceCenter[fi] - cellCentre[cr].
         *		 where cr is the adjoined right cell of face i
         *			   cl is the adjoined left cell of face i.
         */
        virtual void calcFaceCtrToCellCtr() const;

        /**\brief Calculate face weight.*/
        virtual void calcFaceWeight() const;

        /**\brief Calculate node-cell weight.*/
        virtual void calcNodeCellWeight() const;

        /**\brief Calculate cell centres and volumes.*/
        virtual void calcCellCentreAndVol() const;
        void makeCellCentreAndVol(const pointField &p, const vectorArray &fCtrs,
                                  const vectorArray &fAreas, vectorArray &cellCtrs,
                                  realArray &cellVols) const;
        void makeGhostCellCentreAndVol(const vectorArray &fCtrs, const vectorArray &fAreas,
                                       vectorArray &cellCtrs, realArray &cellVols) const;

        /*
         *\brief Vector pointing from the left side cell of face i to the right cell.
         * i.e. adjoinCellCtr[i] = cellCentre[cr] - cellCentre[cl]
         *		where cr is the adjoined right cell of face i
         *			  cl is the adjoined left cell of face i.
         */
        virtual void calcAdjoinnCellCentre() const;

        inline void setGhostCellsType(const short t);

        inline void calcBashMesh() {
#if defined(HUR_FULL_LOGGER)
            if (report) {
                PLDebug("    Debug: calling calcBashMesh().");
            }
#endif // HUR_FULL_LOGGER
            calcFaceCentreAndArea();
            if (!isGhostCellCreated_) {
                createGhostCell();
            }
            modifyFaceZone();
            calcCellCentreAndVol();
            calcFaceWeight();
#if defined(HUR_FULL_LOGGER)
            if (report) {
                PLDebug("    Debug: calling calcBashMesh(). Done.");
            }
#endif // HUR_FULL_LOGGER
        }

        virtual void calcCellAspectRatio() const;

        virtual void calcCellDeltaMax() const;

        void setBaseMesh(const integer nPoint, const integer nFace, const integer nInternalFace,
                         const integer nCell, const bool twoLayer);

        /**\brief Calculate cell weight.*/
        virtual void calcCellWeight(const string &key) const;

        inline void setOriginAndAxis(const point &_origin, const point &_axis) noexcept;

    public:
        declareClassNames;

        baseMesh();

        baseMesh(const controller &cont, const integer nPoint, const integer nFace,
                 const integer nInternalFace, const integer nCell, const bool twoLayer);

        baseMesh(const controller &cont);

        baseMesh(const baseMesh &) = delete;
        baseMesh &operator=(const baseMesh &) = delete;

        /**\brief Destructor.*/
        virtual ~baseMesh() noexcept;

        /**\brief Return the reference of the case control config file.*/
        hur_nodiscard inline const controller &cont() const noexcept;

        /*!\brief Check wether all cells are hexahedral in 3D or quadrilateral in 2D.*/
        void setareAllHexOrQuadri();

        /**\brief Are all cells hexahedral in 3D or quadrilateral in 2D.*/
        hur_nodiscard inline bool areAllHexOrQuadri() const noexcept;

        hur_nodiscard inline short secondNeighbourCellsSize() const noexcept;

        inline void setSecondNeighbourCellsSize(const short s);

        hur_nodiscard inline bool isGhostCellCreated() const noexcept;

        inline void setIsGhostCellCreated(const bool) noexcept;

        /**
         * \brief  The bcType of ghost cell.
         *        -# 0 for no ghost cells;
         *        -# 1 for one layer of ghost cell layers;
         *        -# 2 for two layers of ghost cell layers.
         *        -# 3 for three layers of ghost cell layers, and so on.
         */
        hur_nodiscard inline short ghostCellsType() const noexcept;

        virtual void createGhostCell() = 0;

        /*\brief Can or cannot form two or more ghost cell layers*/
        hur_nodiscard inline bool canTwoLayer() const noexcept { return twoLayer_; }

        /**\brief clear mesh data.*/
        void clearMesh() noexcept;

        /*!\brief Modify faceZone bcType using controller.*/
        void modifyFaceZone();
    };

} // namespace OpenHurricane

#include "baseMesh.inl"