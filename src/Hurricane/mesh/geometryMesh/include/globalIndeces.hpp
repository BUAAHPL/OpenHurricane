/*!
 * \file globalIndeces.hpp
 * \brief Headers of the global indeces.
 *        The subroutines and functions are in the <i>globalIndeces.cpp</i> file.
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
#include "ArrayInclude.hpp"
#include "HurMPI.hpp"
#include "Lists.hpp"
#include "meshElements.hpp"
#include "zoneMesh.hpp"

namespace OpenHurricane {

    // Forward declaration of classes
    class geometryMesh;

    class globalPointIndex

    {
    private:
        const geometryMesh &mesh_;

        integerListList zoneOffset_;

        integerListList processOffset_;

        integerListList pSize_;
        //globalIndex globalIndeces_;

        // (3,10) (5,20)
        std::map<integer, integerVectorList> pointPairMap_;

        std::unordered_map<integer, integer> pointCountMap_;

        /*!\brief Disallow null constructor.*/
        globalPointIndex() = delete;

        /*!\brief set local size for offsets.*/
        void setLocalSize(const integer nLocalPoints);

        /**\bug The number of nonSharedFaceSize used for faceIndexOffSet is
         * not correct.*/
        integer check(const integer nLocalPoints) const;

        void setOffset();

    public:
        /*!\brief Construct from components.*/
        globalPointIndex(const geometryMesh &gMesh,
                         const std::map<integer, integerVectorList> &pointPairMaps);

        /*!\brief Construct from components.*/
        globalPointIndex(const geometryMesh &gMesh,
                         const std::map<integer, integerVectorList> &pointPairMaps,
                         const integerList &pointOffSet);

        /*!\brief Construct from components.*/
        globalPointIndex(const geometryMesh &gMesh);

        /*!\brief Destructor.*/
        ~globalPointIndex() noexcept {}

        // Member Functions

        /*!\brief Clear*/
        inline void clear() noexcept;

        /*!\brief Return mesh reference*/
        hur_nodiscard inline const geometryMesh &mesh() const noexcept;

        void setPointPairMap(const std::map<integer, integerVectorList> &pM);

        /*!\brief Return true if the local point is shared with other process.*/
        hur_nodiscard inline bool isSharedPoint(const integer i) const;

        /*!\brief Return true if the local point is shared with other process.*/
        hur_nodiscard inline bool isSharedPoint(const integer i, bool &isMinId) const;

        /*!\brief Return true if the local point is shared with other process.*/
        hur_nodiscard inline bool isSharedPoint(const integer i, bool &isMinId, int &minId) const;

        /*!\brief Return true if the local point is shared with other process.*/
        inline bool isSharedPoint(const integer i, bool &isMinId, int &minId,
                                  integer &minLocalI) const;

        /*!\brief Return true if the local point is shared with other process*/
        hur_nodiscard inline bool isSharedPoint(const integer i,
                                                integerVectorList &pairsList) const;

        /*!\brief Return true if the local point is shared with other process.*/
        hur_nodiscard inline bool isSharedPoint(const integer i, bool &isMinId,
                                                integerVectorList &pairsList) const;

        /*!\brief Return true if the local point is shared with other process.*/
        hur_nodiscard inline bool isSharedPoint(const integer i, bool &isMinId, int &minId,
                                                integerVectorList &pairsList) const;

        hur_nodiscard inline integer toGlobalIndex(const integer zoneId,
                                                   const integer localI) const;

        hur_nodiscard inline integer toGlobalIndex(const integer localI) const;

        hur_nodiscard inline integer sharedPointSize() const noexcept;

        hur_nodiscard integer whichZone(const integer pi) const;

        /*!\brief Map datatype.*/
        using pointPairMapType = std::map<integer, integerVectorList>;

        void getNodeList(const integer zoneId, pointField &p) const;
    };

    class globalCellIndex {
    private:
        const geometryMesh &mesh_;

        integerListList zoneOffset_;

        integerListList processOffset_;

        void setOffset();

        void checkIndex(const integer zoneId, const integer index) const;

    public:
        /*!\brief Construct from components.*/
        globalCellIndex(const geometryMesh &gMesh);

        inline ~globalCellIndex() noexcept {}

        hur_nodiscard inline integer toGlobalIndex(const integer zoneId, const integer i) const {
            checkIndex(zoneId, i);
            return toGlobalIndex(zoneId, i, HurMPIBase::getProcRank());
        }

        hur_nodiscard inline integer toGlobalIndex(const integer zoneId, const integer i,
                                                   const integer processId) const {
            return zoneOffset_[zoneId][processId] + processOffset_[zoneId][processId] + i;
        }

        hur_nodiscard inline integer toGlobalIndex(const integer i) const {
            auto zoneId = whichZone(i);
            checkIndex(zoneId, i);
            return toGlobalIndex(zoneId, i, HurMPIBase::getProcRank());
        }

        hur_nodiscard integer whichZone(const integer ci) const;
    };

    class globalFaceIndex {
    private:
        const geometryMesh &mesh_;

        integerArrayArray zoneOffset_;

        integerArrayArray processOffset_;

        integerArrayArray cutId_;
        integerArray cutLinkId_;

        /**\brief
         * - first = local face index of currunt procI,
         * - second[0] = remote procI,
         * - second[1] = local face export-id in remote procI,
         * - second[2] = local face id in remote procI,
         * - second[3] = local face's left cell id in remote procI.
         */
        std::map<integer, integerArray> facePairMap_;

        /**\brief Disallow null constructor.*/
        inline globalFaceIndex() = delete;

        /**\brief set local size for offsets.*/
        void setLocalSize(const integer nLocalFace);

        /**\bug The number of nonSharedFaceSize used for faceIndexOffSet is
         * not correct.*/
        integer check(const integer nLocalFace) const;

        void setOffset();

        void getCutZonesIdArrays();

    public:
        /**\brief Construct from components.*/
        globalFaceIndex(const geometryMesh &gMesh,
                        const std::map<integer, integerArray> &facePairMaps);

        /**\brief Construct from compoents.*/
        globalFaceIndex(const geometryMesh &gMesh,
                        const std::map<integer, integerArray> &facePairMaps,
                        const integerArray &faceOffSet);

        /**\brief Construct from components.*/
        globalFaceIndex(const geometryMesh &gMesh);

        /**\brief Destructor.*/
        ~globalFaceIndex() noexcept;

        /**\brief Clear .*/
        inline void clear() noexcept;

        /**\brief Return mesh reference.*/
        hur_nodiscard inline const geometryMesh &mesh() const noexcept;

        void setfacePairMap(const std::map<integer, integerArray> &fM);

        void updatePairMap(const globalCellIndex &gci);

        /**\brief Map datatype.*/
        using facePairMapType = std::map<integer, integerArray>;

        /**\brief Return true if the local face is shared with other process.*/
        hur_nodiscard inline bool isSharedFace(const integer i) const;

        /**\brief Return true if the local face is shared with other process.*/
        hur_nodiscard bool isSharedFace(const integer i, bool &minId) const;

        /**\brief Return true if the local face is shared with other process.*/
        bool isSharedFace(const integer i, bool &isMinId, int &minId) const;

        /**\brief Return true if the local face is shared with other process*/
        hur_nodiscard inline bool isSharedFace(const integer i, integerArray &pairs) const;

        hur_nodiscard integer toGlobalIndex(const integer zoneId, const integer i) const;

        hur_nodiscard integer toGlobalIndex(const integer zoneId, const integer i,
                                            const integer processId) const;

        hur_nodiscard inline integer sharedFaceSize() const noexcept;

        hur_nodiscard inline const facePairMapType &facePairMap() const noexcept;

        void gatherFaceZone(const integer zonei, faceZone &gfz, integerArrayArray &faceConnect,
                            const globalCellIndex &gci, const globalPointIndex &gpi,
                            const faceList &fs, integer &_offset) const;
    };
        
    class globalFaceZoneIndex {
    private:
        /*!\brief Private data.*/

        const geometryMesh &mesh_;

        /**\brief The index of the face zone.*/
        integer id_;

        /**\brief The total points of this face zone.*/
        integer totalNodes_;

        /**\brief The total faces of this face zone.*/
        integer totalFaces_;

        /*
         *\brief The index of points of the face zone to write.
         * The shared point should be written for once.
         */
        integerList writeNodeList_;

        /**\brief The key is the index of the point, and the value is the index of the written nodes.*/
        std::unordered_map<integer, integer> antiOrderMap_;

        /**\brief Only hold in master node.*/
        mutable pointField facePoints_;

        integerList faceDispls_;
        integerList faceRecvcnt_;

        /**\brief Get the node list.*/
        void getNodeList();

        /**\brief Get the node list.*/
        void getNodeListForSerialRun();

        /**\brief Get the points field of the face zone.*/
        void getPoints();

    public:

        /*!\brief Construct from mesh and the index of the face zone.*/
        globalFaceZoneIndex(const geometryMesh &mesh, const integer id);

        /*!\brief Destructor*/
        ~globalFaceZoneIndex() noexcept;

        /**\brief The index of the face zone.*/
        hur_nodiscard inline integer id() const noexcept;

        /**\brief The total points of this face zone.*/
        hur_nodiscard inline integer totalNodes() const noexcept;

        /**\brief The total faces of this face zone.*/
        hur_nodiscard inline integer totalFaces() const noexcept;

        /*!\brief The face zone.*/
        hur_nodiscard const faceZone &fZ() const;

        /*
         *\brief The index of points of the face zone to write.
         * The shared point should be written for once.
         */
        hur_nodiscard inline const integerList &writeNodeList() const noexcept;

        /**\brief The key is the index of the point, and the value is the index of the written nodes.*/
        hur_nodiscard inline const std::unordered_map<integer, integer> &
        antiOrderMap() const noexcept;

        /**\brief Only hold in master node.*/
        hur_nodiscard inline const pointField &facePoints() const noexcept;

        inline void bcastFacePoints() const;

        hur_nodiscard inline const integerList &faceDispls() const noexcept;
        hur_nodiscard inline const integerList &faceRecvcnt() const noexcept;
    };
} //  namespace OpenHurricane

#include "globalIndeces.inl"