/*!
 * \file geometryMesh.hpp
 * \brief Headers of the geometry mesh.
 *        The subroutines and functions are in the <i>geometryMesh.cpp</i> file.
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

#include "CRSMatrixAddressing.hpp"
#include "baseMesh.hpp"
#include "dataStructure.hpp"
#include "globalIndeces.hpp"
#include "processTopologies.hpp"
#include "registerTable.hpp"
#include "tecplotFormat.hpp"
#include "zoneMesh.hpp"

namespace OpenHurricane {

    // Forward declaration of classes
    class globalMesh;
    class originMeshRead;
    class decomposingMeshByMetis;
    class hdf5O;
    class hdf5I;

    class geometryMeshProcs {
    protected:
        integerArray partOfProcs_;
        integerArrayArray cellOfProc_;
        integerArray perm_;

        integerArray iperm_;

    public:
        inline geometryMeshProcs() : partOfProcs_(), cellOfProc_(), perm_(), iperm_() {}
        geometryMeshProcs(const geometryMeshProcs &) = delete;
        geometryMeshProcs &operator=(const geometryMeshProcs &) = delete;

        inline ~geometryMeshProcs() noexcept {}

        hur_nodiscard inline const integerArray &partOfProcs() const noexcept {
            return partOfProcs_;
        }

        hur_nodiscard inline const integerArrayArray &cellOfProc() const noexcept {
            return cellOfProc_;
        }

        hur_nodiscard inline const integerArray &perm() const noexcept { return perm_; }

        hur_nodiscard inline const integerArray &iperm() const noexcept { return iperm_; }
    };

    class geometryMeshCore : public geometryMeshProcs {
    public:
        using originIdMapType = std::map<integer, integer>;

    protected:
        /*!\brief Origin cell index of .h5 file mesh.*/
        integerArrayArray originCellIndex_;

        /*\!brief Map of origin mesh index and relay file index.
         * Note: key is origin mesh index, value is relay file index.
         */
        originIdMapType indexMap_;

        /*!\brief The number of points of the full mesh.*/
        integer allMeshNodeNumber_;

        /*!\brief The number of faces of the full mesh.*/
        integer allMeshFaceNumber_;

        /*!\brief The number of cells of the full mesh.*/
        integer allMeshCellNumber_;

        mutable uniquePtr<integer> totalCutFacesPtr_;

        mutable uniquePtr<integer> totalPerFacesPtr_;

    public:
        inline geometryMeshCore()
            : geometryMeshProcs(), originCellIndex_(), indexMap_(), allMeshNodeNumber_(0),
              allMeshFaceNumber_(0), allMeshCellNumber_(0), totalCutFacesPtr_(nullptr),
              totalPerFacesPtr_(nullptr) {}
        geometryMeshCore(const geometryMeshCore &) = delete;
        geometryMeshCore &operator=(const geometryMeshCore &) = delete;

        inline ~geometryMeshCore() noexcept {}

        /*!\brief Origin cell index of .h5 file mesh.*/
        hur_nodiscard inline const integerArrayArray &originCellIndex() const noexcept {
            return originCellIndex_;
        }

        /*\!brief Map of origin fluent index and relay file index.*/
        hur_nodiscard inline originIdMapType &indexMap() noexcept { return indexMap_; }

        /*\!brief Map of origin fluent index and relay file index.*/
        hur_nodiscard inline const originIdMapType &indexMap() const noexcept { return indexMap_; }

        hur_nodiscard inline integer allMeshNodeNumber() const noexcept {
            return allMeshNodeNumber_;
        }

        hur_nodiscard inline integer allMeshFaceNumber() const noexcept {
            return allMeshFaceNumber_;
        }

        hur_nodiscard inline integer allMeshCellNumber() const noexcept {
            return allMeshCellNumber_;
            ;
        }
    };

    /**
     * \brief The class of geometry mesh.
     */
    class geometryMesh : public baseMesh, public registerTable, public geometryMeshCore {
    public:
        using Mesh = geometryMesh;

    private:
        void getCellOfProc();

        void transferFromDCM(const controller &cont, decomposingMeshByMetis &dcM);

        /*!\brief Get scale factor from mesh unit.
         * \param[in] meshUnit - The unit of mesh.
         * \return Return the scale factor. (1.0 for m).
         */
        hur_nodiscard real meshScaleFactor(const string &meshUnit) const;

    private:
        void calcNumberOfCutFaces() const;

        void calcNumberOfPerFaces() const;

        /*!\brief Parallel info.*/
        mutable uniquePtr<globalMesh> globalMeshInfoPtr_;

        /**\brief The global face zone information.*/
        mutable sharedPtrList<globalFaceZoneIndex> globalFaceZoneL_;

        /*!\brief Periodic info*/
        std::map<integer, integerArray> perFaceMaps_;

        /*!
         *\brief Number of pairs of process to process.
         * That is the size of the cut zone list from decomposing mesh.
         */
        integer pairProcessSize_;

        /*!
         *\brief Number of pairs of process to process.
         * That is the size of the periodic zone list from decomposing mesh.
         */
        integer pairPeriodicProcessSize_;

        integer periodicPairSize_;

        /** \brief Cell Orthogonality. */
        mutable uniquePtr<realArray> cellOrthoPtr_;

        void calcCellOrtho(const realArray &faceOrtho, realArray &cellOrtho) const;
        void makeCellOrtho(const realArray &faceOrtho) const;
        void makeCellOrtho() const;

        // For Face skewness correction.

        /**
         * \brief The intersection point of the line l_cl_cr and the face.
         *  l_cl_cr = cellCentre[cl]-cellCentre[cr]
         */
        mutable uniquePtr<vectorArray> faceCellIntersectionPtr_;
        void makeFaceCellIntersectionPoint() const;

        /**
         * \brief The vector pointing from the intersection point of the line l_cl_cr and the face to the face centre.
         */
        mutable uniquePtr<vectorArray> fIntSectFCentrePtr_;
        void makeFIntSectFCentre() const;

        /** \brief Face weight for skewness correction. */
        mutable uniquePtr<realArray> faceSkewnessWeightPtr_;
        void makeFaceSkewnessWeight(const real minSkewness) const;

        mutable uniquePtr<integerArray> skewFacePtr_;
        void makeSkewFace(const real minSkewness) const;

        mutable uniquePtr<boolList> isSkewFacePtr_;
        void makeIsSkewFace(const real minSkewness) const;

        // LDU addressing for LDU matrix
        mutable uniquePtr<integerArray> LDUFaceMapPtr_;

        void calcLDUFaceMap() const;

    private:
        mutable uniquePtr<CRSMatrixAddressing> CRSAdrrPtr_;
        void calcCRSAdrr() const;

    public:
        hur_nodiscard inline const CRSMatrixAddressing &CRSAdrr() const {
            if (!CRSAdrrPtr_) {
                calcCRSAdrr();
            }
            return *CRSAdrrPtr_;
        }

    private:
        /*!\brief Create first layer cut zone.*/
        void createFirstLayerCutZone(const std::map<integer, integerArray> &_faceMap);

        void createFirstLayerPerZone(const std::map<integer, integerArray> &_faceMap);

        /*!\brief Create second or more layers ghost cells.*/
        void createSecondGhostCellLayer(const short nLayer,
                                        const std::map<integer, integerArray> &_faceMap,
                                        const std::map<integer, integerArray> &_facePeriMap,
                                        cutZoneList &cZ, perZoneList &pZ);

        // process topology

    private:
        mutable uniquePtr<processCutTopology> processCutPtr_;
        mutable uniquePtr<processPerTopology> processPerPtr_;

        void makeProcessCutTopology() const;
        void makeProcessPerTopology() const;

    public:
        hur_nodiscard inline const processCutTopology &processCut() const {
            if (!processCutPtr_) {
                makeProcessCutTopology();
            }
            return *processCutPtr_;
        }
        hur_nodiscard inline const processPerTopology &processPer() const {
            if (!processPerPtr_) {
                makeProcessPerTopology();
            }
            return *processPerPtr_;
        }

    protected:
        hur_nodiscard fileName getMeshFileName(const controller &cont) const;

        void readMesh(const controller &cont);

        void setCoupledWall(const controller &cont, originMeshRead &gmRead);

        /*!\brief Scale mesh with a scalar factor.*/
        void scaleMesh(const real scaleFactor);

        /*!\brief Scale mesh with different factor on different direction.*/
        void scaleMesh(const vector &scaleFactor);

        //void setWriteControl();

        void writeASCIIToBinary(const char *str, fileOsstream &fos);

        integer getPeriodicPairSize(const perZoneList &pZL, integerArray &pi) const;

        void getPeriodicPairList(const perZoneList &pZL, periodicPairList &ppl) const;
        void getPeriodicPair(const perZone &pZL, periodicPair &ppl) const;

    public:
        /**\brief Return a null geometry mesh.*/
        hur_nodiscard inline static const geometryMesh &nullObject();

        geometryMesh(const object &ob);

        geometryMesh(object &&ob);

        /**\brief Construct from components.*/
        geometryMesh(const object &ob, const controller &cpnt, pointField &meshPoints,
                     faceList &meshFaces, const integer nInterFaces, cellList &meshCells,
                     pointZoneList &pointZones_, faceZoneList &faceZones_, cutZoneList &cutZones_,
                     perZoneList &perZones_, cellZoneList &cellZones_, const bool twoLayer);

        geometryMesh(object &&ob, const controller &cont, pointField &meshPoints,
                     faceList &meshFaces, const integer nInterFaces, cellList &meshCells,
                     pointZoneList &pointZones_, faceZoneList &faceZones_, cutZoneList &cutZones_,
                     perZoneList &perZones_, cellZoneList &cellZones_, const bool twoLayer);

        /*!\brief Construct from decomposedMesh*/
        geometryMesh(const object &ob, const controller &cont);

        /*!\brief Construct from decomposedMesh*/
        geometryMesh(object &&ob, const controller &cont);

        /*!\brief Construct from decomposedMesh*/
        geometryMesh(object &&ob, const controller &cont, const std::string &meshStr);


        geometryMesh(const geometryMesh &) = delete;
        geometryMesh &operator=(const geometryMesh &) = delete;

        /**\brief Destructor.*/
        virtual ~geometryMesh() noexcept;

        /*!
         * \brief Return size of the internal field
         *        (for geometry, is the size of internal cells).
         */
        hur_nodiscard inline integer internalArraySize() const;

        /*!\brief Return the total cells of this mesh.*/
        hur_nodiscard inline integer size() const noexcept;

        /*!\brief Is the ith face an interior face.*/
        hur_nodiscard inline bool isInteriorFace(const integer i) const noexcept;

        /*!\brief Is the ith face a boundary face.*/
        hur_nodiscard inline bool isBoundaryFace(const integer i) const noexcept;

        /*!\brief It the ith cell an internal cell.*/
        hur_nodiscard inline bool isInternalCell(const integer i) const noexcept;

        /*!\brief It the ith cell an dummy cell.*/
        hur_nodiscard inline bool isDummyCell(const integer i) const noexcept;

        /**\brief Return the total number of the cut face.*/
        hur_nodiscard inline integer totalCutFaces() const;

        /**\brief Return the total number of the periodic face.*/
        hur_nodiscard inline integer totalPerFaces() const;

        /*!\brief Create ghost cells for boundary.*/
        virtual void createGhostCell();

        /**\brief Return the number of face zones except cut zones*/
        hur_nodiscard integer nFaceZoneExceptCutZone() const;

        // Check mesh

        void checkMeshFaceSkewness(real &maxSkew) const;

        void checkMeshOrthogonality(real &minOrtho) const;

        void checkMeshVolRatio(real &maxRatio) const;

        void checkMeshAspectRatio(real &maxRatio) const;

        void checkAndRportMeshQuality() const;

        // Global mesh

        hur_nodiscard const globalMesh &globalMeshInfo() const;

        hur_nodiscard const globalFaceZoneIndex &globalFaceZoneInfo(const integer fzid) const;

        /** \brief Cell Orthogonality. */
        hur_nodiscard const realArray &cellOrtho() const;

        /**
         * \brief The intersection point of the line l_cl_cr and the face.
         *  l_cl_cr = cellCentre[cl]-cellCentre[cr]
         */
        hur_nodiscard const vectorArray &faceCellIntersectionPoint() const;

        /**
         * \brief The vector pointing from the intersection point of the line l_cl_cr and the face to the face centre.
         */
        hur_nodiscard const vectorArray &fIntSectFCentre() const;

        /** \brief Face weight for skewness correction. */
        hur_nodiscard const realArray &faceSkewnessWeight(const real minSkewness = 0.1) const;

        hur_nodiscard const integerArray &skewFace(const real minSkewness = 0.1) const;

        hur_nodiscard const boolList &isSkewFace(const real minSkewness = 0.1) const;

        /**\brief Is mesh moving (It is not used).*/
        hur_nodiscard inline bool moving() const noexcept { return false; }

        /*!\brief Clear the mesh.*/
        void clearGeoMesh() noexcept;

        /**
         * \brief Face map from geometry mesh to LDU matrix.
         */
        hur_nodiscard const integerArray &LDUFaceMap() const;

        // Write

        /*!\brief Write mesh to tecplot.*/
        void writeMeshToTecplot(fileOsstream &fos) const;

        /*!\brief Write mesh of face zone with index: fzid to tecplot.*/
        void writeMeshToTecplot(fileOsstream &fos, const integer fzid) const;

        /*!
         * \brief Write tecplot file header.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         */
        void writeTecplotHeader(fileOsstream &fos,
                                const short fileType = tecplotFormat::FULL) const;

        /*!
         * \brief Write tecplot file header.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         */
        void writeTecplotHeader(fileOsstream &fos, const stringList &outVarName,
                                const short fileType = tecplotFormat::FULL) const;

        /*!
         * \brief Write tecplot file header for face zone with index: fzid.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         */
        void writeTecplotHeader(fileOsstream &fos, const integer fzid,
                                const short fileType = tecplotFormat::FULL) const;

        /*!
         * \brief Write tecplot file header for face zone with index: fzid.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         */
        void writeTecplotHeader(fileOsstream &fos, const integer fzid, const stringList &outVarName,
                                const short fileType = tecplotFormat::FULL) const;

        /*!
         * \brief Write zone header of tecplot file.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         * \param[in] faceConnect - True if it is a parallel run for cell-based result output.
         */
        void writeTecplotZoneHeader(fileOsstream &fos, const short fileType = tecplotFormat::FULL,
                                    bool faceConnect = false) const;

        /*!
         * \brief Write zone header of tecplot file.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         * \param[in] faceConnect - True if it is a parallel run for cell-based result output.
         */
        void writeTecplotZoneHeader(fileOsstream &fos, const stringList &outVarName,
                                    const short fileType = tecplotFormat::FULL,
                                    bool faceConnect = false) const;

        /*!
         * \brief Write zone header of tecplot file for face zone with index: fzid.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         */
        void writeTecplotZoneHeaderMaster(fileOsstream &fos, const integer fzid,
                                          const short fileType = tecplotFormat::FULL) const;

        /*!
         * \brief Write zone header of tecplot file for face zone with index: fzid.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         * \param[in] faceConnect - True if it is a parallel run for cell-based result output.
         */
        void writeTecplotZoneHeaderMaster(fileOsstream &fos, const integer fzid,
                                          const stringList &outVarName,
                                          const short fileType = tecplotFormat::FULL) const;

        /*!\brief Write the points of mesh to tecplot file.*/
        void writePointToTecplot(fileOsstream &fos) const;

        /*!\brief Write the points of face zone with index: fzid to tecplot file.*/
        void writePointToTecplot(fileOsstream &fos, const integer fzid) const;

        /*!\brief Write cell connectivity to tecplot file.*/
        void writeCellToTecplot(fileOsstream &fos) const;

        /*!\brief Write face connectivity of face zone with index: fzid to tecplot file.*/
        void writeFaceToTecplot(fileOsstream &fos, const integer fzid) const;

        /*!\brief Write face connectivity to tecplot for parallel running case.*/
        void writeFaceConnectTecplot(fileOsstream &fos) const;

        /*!
         * \brief Write result to tecplot file.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         */
        void writeOutputToTecplot(fileOsstream &fos,
                                  const short fileType = tecplotFormat::FULL) const;

        /*!
         * \brief Write result to tecplot file.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         */
        void writeOutputToTecplot(fileOsstream &fos, const stringList &outVarName,
                                  const short fileType = tecplotFormat::FULL) const;

        /*!
         * \brief Write result to tecplot file.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         */
        void writeOutputToTecplot(fileOsstream &fos, const integer fzid,
                                  const short fileType = tecplotFormat::FULL) const;

        /*!
         * \brief Write result to tecplot file.
         * \param[in] fileType
         * - FULL: mesh and solution;
         * - GRID: only mesh;
         * - SOLUTION: only solution.
         */
        void writeOutputToTecplot(fileOsstream &fos, const integer fzid,
                                  const stringList &outVarName,
                                  const short fileType = tecplotFormat::FULL) const;

        void writeMesh(hdf5O &fos, const string &gridGroupName) const;

    private:
        void writePoint(hdf5O &fos, const string &gridGroupName) const;

        void writeCell(hdf5O &fos, const string &gridGroupName) const;

        void writeFace(hdf5O &fos, const string &gridGroupName) const;

        void writePeriodicPairList(hdf5O &fos, const string &gridGroupName) const;

        void writePointZoneAttribute(hdf5O &fos, const pointZone &pz, const pointField &p,
                                     const string &gridGroupName, const string &dataName) const;

        void writeCellZoneAttribute(hdf5O &fos, const cellZone &cz, const string &gridGroupName,
                                    const string &dataName) const;

        void writeFaceZoneAttribute(hdf5O &fos, const faceZone &fz, const string &gridGroupName,
                                    const string &dataName) const;

        void writePeriodicPairAttribute(hdf5O &fos, const periodicPair &fz,
                                        const string &gridGroupName, const string &dataName) const;

        hur_nodiscard integerArrayArray getCutZonesIdArrays() const;

        hur_nodiscard pointField gatherAllPoints() const;

    public:
        // Find points

        /**
         * \brief Find the cell that contains the given point.
         * \param[in] p - The given point
         * \return The cell index (return -1 if not found)
         * \retval An integer value
         */
        hur_nodiscard integer findPointCell(const point &p) const;

        hur_nodiscard bool pointInCell(const point &p, const integer celli) const;

        hur_nodiscard vectorArray allCellCentre() const;

        bool getInterpolationSource(const hdf5I &fos);

    protected:
        /*!\brief Cell centre of source interpolation mesh.*/
        mutable uniquePtr<vectorArray> sorCellCentrePtr_;

        /*!\brief Target neighbour cell of source mesh.*/
        mutable uniquePtr<integerArrayArray> sorKNN_;

        mutable uniquePtr<integerArrayArray> tarOfProc_;

    public:
        hur_nodiscard inline const integerArrayArray &tarOfProc() const noexcept;

        hur_nodiscard inline const vectorArray &sorCellCentre() const noexcept;

        hur_nodiscard inline const integerArrayArray &sorKNN() const noexcept;

    protected:
        /**
         * \brief Load weights of each cells.
         */
        mutable uniquePtr<integerArray> cellLoadWeightsPtr_;

        void writeCellLoadWeights(hdf5O &fos, const string &gridGroupName) const;

    public:
        /**
         * \brief Load weights of each cells.
         */
        hur_nodiscard inline const integerArray &cellLoadWeights() const noexcept;

        /**
         * \brief Load weights of each cells.
         */
        hur_nodiscard inline integerArray &cellLoadWeights() noexcept;

        hur_nodiscard inline bool hasCellLoadWeights() const noexcept;
    };

    class meshCheckers {
    public:
        static void minMaxFaceArea(const vectorArray &faceArea, const faceZoneList &fZL,
                                   real &minArea, real &maxArea);

        static void minMaxCellVol(const realArray &cellVol, const cellZoneList &cZL, real &minVol,
                                  real &maxVol, real &totalVol);

        static void minMaxCellCentre(const pointField &cellCntr, const cellZoneList &cZL,
                                     vector &mincCComponents, vector &maxcCComponents);

        static void minMaxPointDomain(const pointField &point, const integer nPoints,
                                      vector &minPointCompo, vector &maxPointCompo);

        /*!\brief To generate the face skewness field.*/
        static hur_nodiscard realArray faceOrthogonality(const geometryMesh &mesh,
                                                         const vectorArray &faceAreas,
                                                         const vectorArray &faceCentres,
                                                         const vectorArray &cellCentres);

        /*!\brief To generate the face skewness field.*/
        static hur_nodiscard realArray faceSkewness(const geometryMesh &mesh,
                                                    const pointField &points,
                                                    const vectorArray &faceCentres,
                                                    const vectorArray &faceAreas,
                                                    const vectorArray &cellCentres);

        /*!\brief To generate the volume ratio field.*/
        static hur_nodiscard realArray volRatio(const geometryMesh &mesh, const realArray &volume);
    };

} //  namespace OpenHurricane

#include "geometryMesh.inl"