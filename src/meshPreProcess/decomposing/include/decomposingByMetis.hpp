/*!
 * \file decomposingByMetis.hpp
 * \brief Headers of decomposing mesh by metis.
 *        The subroutines and functions are in the <i>decomposingByMetis.cpp</i> file.
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

#include "Lists.hpp"
#include "geometryMesh.hpp"
#include "loadBalancingWeights.hpp"
#include "logFile.hpp"
#include "originMeshRead.hpp"

namespace OpenHurricane {
    class decomposingMeshByMetis {
    private:
        const originMeshRead &originMesh_;

        vector origin_;
        vector axis_;

        integerArray partOfProcs_;

    public:
        inline const integerArray &partOfProcs() const noexcept;
        inline integerArray &partOfProcs();

    private:
        integerArray xadj_;

        integerArray adjncy_;

        /*!\brief The number of partitions of the full mesh.*/
        integer nparts_;

        // Information of originMesh decomposed to current processor

        /*!\brief Information of cell in current processor.*/
        cellList cells_;

        /*!\brief Information of face in current processor.*/
        faceList faces_;

        /*!\brief Information of node in current processor.*/
        pointField nodes_;

        /*!\brief Information of cellZone in current processor.*/
        cellZoneList cellZones_;

        /*!\brief Information of faceZone in current processor.*/
        faceZoneList faceZones_;

        /*!\brief Information of pointZone in current processor.*/
        pointZoneList pointZones_;

        /*!\brief Information of cutZone in current processor
         * Only assigned in master processor.*/
        cutZoneList cutZones_;

        /*!\brief Information of face shared by different processores
         * include interiorwall faces.*/
        cutZoneList faceShareZones_;

        /*!\brief Information of perZone in current processor.*/
        perZoneList perZones_;

        integer periodicPairSize_;

        /*!\brief Information of facePairMaps in current processor.*/
        /*!\brief integerArray contains 4 elements: recvP, exportIndex, localIndex, localCellIndex.*/
        //std::map<integer, integerVector2D> facePairMaps_;
        //std::map<integer, integerVector> facePairMaps_;
        std::map<integer, integerArray> facePairMaps_;

        /*!\brief Information of pointPairMaps in current processor.*/
        std::map<integer, integerVectorList> pointPairMaps_;

        /**\brief Information of periodic facePairMaps in current processor*/
        std::map<integer, integerArray> perPairMaps_;

        /*!\brief Information of cell offset.*/
        integerArray cellOffSet_;

        /*!\brief Information of face offset.*/
        integerArray faceOffSet_;

        /*!\brief Information of point offset.*/
        integerArray pointOffSet_;

        /*!\brief Information of cell local index in master processor.*/
        integerArray color_;

        /*!\brief Information of whether originMesh_ is structed or not.*/
        bool unStructed_;

        /*!\brief Information of whether originMesh_ has structed cellZone or not.*/
        bool hasStructed_;

        /*!\brief Form two ghost layers or not.*/
        bool twoGhostLayer_;

        /*!\brief Disallow null constructor.*/
        decomposingMeshByMetis() = delete;

        /*!\brief Disallow construct as copy.*/
        decomposingMeshByMetis(const decomposingMeshByMetis &) = delete;

        /*!\brief Disallow default bitwise assignment.*/
        void operator=(const decomposingMeshByMetis &) = delete;

        void computeArray();

        void clear() noexcept;

        void cellToProcessor();

        void faceToProcessor(integerArrayArray &);

        void PointToProcessor(const integerArrayArray &);

        void PointToProcessorFix(const integerArrayArray &faceOfProc);

        //void perFaceToProcessor(integerArrayArray&, integerArrayArrayArray&, integerArray&);
        void perFaceToProcessor(integerArray &, integerArrayArray &);

        template <class T> void bcastZone(List<T> &);

        void cellAlloc(integerArrayArray &);

        void cellDistribute(integerArrayArray &, /*integerArray&,*/ integer *, const bool,
                            const integer);

        void cellScatter(integerArrayArray & /*, integerArray&*/);

        void cellMember(const integerArrayArray &, integerArray &, integerArray &, integer *,
                        integer *, integer *);

        void cellReconstr(const integerArray &, const integerArray &, const integerArray &,
                          const integer);

        void perZoneAlloc(integerArrayArray &, integerArray &);

        //void perZoneScatter(integerArrayArray&, integerArrayArray&, integerArrayArrayArray&, integerArray&);
        void perZoneScatter(integerArrayArray &, const integerArrayArray &, integerArray &);

        /*void perFaceDistribute
        (
                const periodicPair&, integerArrayArray&, integerArrayArray&, integerArrayArray&,
                integerArrayArrayArray&, integerArray&, const integerArrayArray&, const bool
        );*/
        void perFaceDistribute(const periodicPair &, integerArray &, integerArrayArray &,
                               integerArrayArray &, const integerArrayArray &, const bool);

        void createPerZone(const periodicPair &, integerArray &, integerArrayArray &,
                           integerArrayArray &, const integerArrayArray &,
                           const integerArrayArray &, perZoneList &newPerZone);

        void faceAlloc(integerArrayArray &, integerArrayArray &);

        void faceScatter(integerArrayArray &, integerArrayArray &, integerArrayArrayArray &,
                         integerArray &);

        void faceDistribute(integerArray &, integerArrayArray &, integerArrayArray &,
                            integerArrayArray &, integerArrayArrayArray &, integerArray &,
                            const bool, const bool, const integer);

        void faceMember(const integer, const integerArrayArray &, integerArrayArrayArray &,
                        integerArray &, List<short> &, integerArray &, integerArray &,
                        integerArray &);

        void faceReconstr(const integerArray &, const integerArray &, const List<short> &,
                          const integer);

        void checkPartition(const integerArrayArrayArray &);

        //void faceShare(const integerArray&);
        void faceShare(const integerArrayArray &, const integerArray &, const integer);

        void pointShare(const integerArrayArray &, integerArray &, const integerArrayArray &,
                        const integer);

        /**
         * \brief To find shared points.
         * \param[in] faceOfProc - The face list of each process.
         * \param[out] sizeLP - The size of points at each process
         * \param[in] mapOP - Map of point. Each global point corresponds to a list mapOP[pointGlobaleIndex],
         *                     of which each element stores integerVectod2D (ip,localId), meaning that the point is shared in ip process with the local index (localId).
         * \param[in] nFace - The number of faces of the whole mesh.
         */
        void pointShareFix(const integerArrayArray &faceOfProc, integerArray &sizeLP,
                           const List<List<Vector2D<integer>>> &mapOP, const integer nFace);

        void share(cutZoneList &, const integerArrayArray &, const integerArray &,
                   const integerArrayArray &, const integer, integerArrayArray &,
                   integerArrayArray &, const bool);

        /**
         * \brief To find shared points.
         * \param[in] faceOfProc - The face list of each process.
         * \param[out] sizeLP - The size of points at each process
         * \param[in] mapOP - Map of point. Each global point corresponds to a list mapOP[pointGlobaleIndex],
         *                     of which each element stores integerVectod2D (ip,localId), meaning that the point is shared in ip process with the local index (localId).
         * \param[in] nFace - The number of faces of the whole mesh.
         */
        void shareFix(cutZoneList &cutZones, const integerArrayArray &faceOfProc,
                      const integerArray &sizeLP, const List<List<Vector2D<integer>>> &mapOP,
                      const integer nFace, integerArrayArray &globalColor, integerArrayArray &cZone,
                      const bool pShare);

        void pointMap(List<List<std::map<integer, integer>>> &, integerArray &,
                      const integerArrayArray &, const integerArrayArray &);

        void pointMapFix(List<List<std::map<integer, integer>>> &mapColor, integerArray &mapSize,
                         const integerArrayArray &faceOfProc,
                         const List<List<Vector2D<integer>>> &mapOP);

        void faceMap(const integer, const integerArrayArray &, const integerArrayArray &,
                     std::map<integer, integerArray> &);

        template <class T> void formCutPerZone(List<T> &);

        void isUnStructedMesh();

        template <class T>
        void isTwoLayer(const integer, const integer, const integer, const integer,
                        const integerArrayArray &, List<T> &);

        /*! Check whether the geometryMesh do can form two ghost layer*/
        bool isTwoLayer();

        void pointScatter(const integerArrayArray &, integerArrayArray &, integerArray &);

        /**
         * \brief To scatter points to other processes.
         * \param[in] faceOfProc - The face list of each process.
         * \param[out] mapOP - Map of point. Each global point corresponds to a list mapOP[pointGlobaleIndex],
         *                     of which each element stores integerVectod2D (ip,localId), meaning that the point is shared in ip process with the local index (localId).
         * \param[out] sizeLP - The size of points at each process
         */
        void pointScatterFix(const integerArrayArray &faceOfProc,
                             List<List<Vector2D<integer>>> &mapOP, integerArray &sizeLP);

        void cellAssign();

        void getPeriodicTransformVector(const periodicPair &pp, point &tranfV) const;

        void getPeriodicRotationalMatrix(const periodicPair &pp, tensor &RM) const;

        void calcRotationMatrix(const vector ra, const real theta, tensor &RM) const;

        void createFaceShareZone(const integerArray &globalColor);

    private:
        /**\brief The lists of the processors used in the origin mesh (the input mesh).*/
        integerArrayArray decomposeList_;
        /**\brief The number of the processors used in the origin mesh (the input mesh).*/
        integer originMeshDecompSize_;

        void decomposeFromOriginMesh();

    public:
        /*!\brief Construct from components.*/
        decomposingMeshByMetis(const originMeshRead &originMeshes, const integer nP,
                               const controller &cont);

        /*!\brief Destructor.*/
        inline ~decomposingMeshByMetis() noexcept;

        inline const originMeshRead &originMesh() const;

        inline integer nProcessors() const;

        inline bool isUnStructed() const;

        inline bool hasStructed() const;

        inline bool twoGhostLayer() const;

        inline const vector &origin() const noexcept;
        inline const vector &axis() const noexcept;

        inline cellList &cells();

        inline faceList &faces();

        inline pointField &points();

        inline cellZoneList &cellZones();

        inline faceZoneList &faceZones();

        inline pointZoneList &pointZones();

        inline perZoneList &perZones();

        inline cutZoneList &cutZones();

        inline integer periodicPairSize() const noexcept;

        inline const std::map<integer, integerArray> &facePairMaps() const;

        inline std::map<integer, integerArray> &facePairMaps();

        inline std::map<integer, integerArray> &perPairMaps();

        inline const std::map<integer, integerVectorList> &pointPairMaps() const;

        inline const integerArray &cellOffSet() const;

        inline const integerArray &faceOffSet() const;

        inline const integerArray &pointOffSet() const;

        integer nInteriorFaces() const;

        // Decompose gobal mesh to local mesh and send

        /*!\brief Mark to which processor a cell belong.*/
        int decomposing();

    protected:
        /*!\brief Mark to which processor a cell belong by given the weights of the cells.*/
        int decomposing(const integerArray &cellWeight);

    public:
        void meshToProcessor();

    protected:
        // Load balancing weights
        loadBalancingWeights *loadWgtPtr_;
    };
} // namespace OpenHurricane

template <class Type> void OpenHurricane::decomposingMeshByMetis::bcastZone(List<Type> &zoneList) {
    integer nzone;
    nzone = zoneList.size();
    HurMPI::bcast(&nzone, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

    integer *index;
    integer *charSize;
    short *types;

    index = new integer[nzone];
    charSize = new integer[nzone];
    types = new short[2 * nzone];

    std::string name;
    integer nameSize = 0;

    if (HurMPI::master()) {
        for (integer i = 0; i < nzone; i++) {
            index[i] = zoneList[i].index();
            feature<Type>::pack(zoneList[i], types + 2 * i);
            name += zoneList[i].name();
            charSize[i] = integer(zoneList[i].name().size());
        }
        nameSize = integer(name.size());
    }

    HurMPI::bcast(&nameSize, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

    char *names = new char[nameSize + 1]();

    if (HurMPI::master()) {
        strcpy(names, name.c_str());
    }
    names[nameSize] = '\0';

    HurMPI::bcast(index, nzone, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
    HurMPI::bcast(names, nameSize + 1, MPI_CHAR, HurMPI::masterNo(), HurMPI::getComm());
    HurMPI::bcast(charSize, nzone, feature<integer>::MPIType, HurMPI::masterNo(),
                  HurMPI::getComm());
    HurMPI::bcast(types, 2 * nzone, feature<short>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

    List<Type> zones(nzone);
    if (!HurMPI::master()) {
        zoneList.resize(nzone);
        feature<Type>::objConstruct(zones, index, charSize, types, names, nzone);
        zoneList.transfer(zones);
    }

    delete[] index;
    delete[] charSize;
    delete[] types;
    delete[] names;
}

template <class Type>
void OpenHurricane::decomposingMeshByMetis::isTwoLayer(const integer ip1, const integer ip2,
                                                       const integer cell1, const integer cell2,
                                                       const integerArrayArray &cpColor,
                                                       List<Type> &cpZones) {
    if (cpZones[cpColor[ip1][ip2]].ghostLayer()) {
        if (isUnStructed()) {
            if (hasStructed()) {
                for (integer zoneI = 0; zoneI < originMesh_.cellZones().size(); ++zoneI) {
                    if ((cell1 >= originMesh_.cellZones()[zoneI].firstIndex() &&
                         cell1 <= originMesh_.cellZones()[zoneI].lastIndex() &&
                         originMesh_.cellZones()[zoneI].shapeType() !=
                             cellShapeType::shapeTypes::hexahedral) ||
                        (cell2 >= originMesh_.cellZones()[zoneI].firstIndex() &&
                         cell2 <= originMesh_.cellZones()[zoneI].lastIndex() &&
                         originMesh_.cellZones()[zoneI].shapeType() !=
                             cellShapeType::shapeTypes::hexahedral)) {
                        cpZones[cpColor[ip1][ip2]].setGhostLayer(false);
                        cpZones[cpColor[ip2][ip1]].setGhostLayer(false);
                    }
                }
            } else {
                cpZones[cpColor[ip1][ip2]].setGhostLayer(false);
                cpZones[cpColor[ip2][ip1]].setGhostLayer(false);
            }
        }
    }
}

template <class Type>
void OpenHurricane::decomposingMeshByMetis::formCutPerZone(List<Type> &cpZone) {
    for (integer zoneI = 0; zoneI < cpZone.size(); ++zoneI) {
        integer sendP = cpZone[zoneI].sendProc();
        if (HurMPI::isThisProc(sendP)) {
            integer recvP = cpZone[zoneI].receivProc();
            integer size = cpZone[zoneI].sor().size();
            integerArray sorF(size);
            integerArray desF(size);
            for (integer i = 0; i < cpZone[zoneI].sor().size(); ++i) {
                sorF[i] = faces_[cpZone[zoneI].sor()[i]].leftCell();
                desF[i] = cpZone[zoneI].des()[i] + cellOffSet_[recvP];
            }
            //if (cpZone[zoneI].ghostLayer())
            if (twoGhostLayer()) {

                integerArray sorS(size);
                integerArray desS(size);
                for (integer i = 0; i < cpZone[zoneI].sor().size(); ++i) {
                    integer cellI = sorF[i];
                    integer faceI = cells_[cellI].findFace(cpZone[zoneI].sor()[i]);

                    std::map<short, short>::const_iterator iter;
                    iter = cellFacePair::HexCellFacePairMap.find(faceI);
                    if (iter == cellFacePair::HexCellFacePairMap.end()) {
                        LFatal("Not find opposite face, maybe cell type error ");
                    }
                    integer oppositeFace = cells_[cellI].facesList()[iter->second];
                    integer oppositeCell =
                        faces_[oppositeFace].leftCell() + faces_[oppositeFace].rightCell() - cellI;
                    sorS[i] = oppositeCell;
                    desS[i] = desF[i] + size;
                }
                sorF.append(sorS);
                desF.append(desS);
            }
            cpZone[zoneI].sor().resize(sorF.size());
            cpZone[zoneI].des().resize(desF.size());
            cpZone[zoneI].sor() = sorF;
            cpZone[zoneI].des() = desF;
        }
    }
}

#include "decomposingByMetis.inl"