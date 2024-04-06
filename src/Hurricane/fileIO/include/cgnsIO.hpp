/*!
 * \file cgnsIO.hpp
 * \brief Headers of input and output cgns file.
 *        The subroutines and functions are in the <i>cgnsIO.cpp</i> file.
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
#include <iostream>
#ifdef USES_CGNS

#ifdef WIN32
#include "cgnslib.h"
#include "cgnstypes.h"
#else
#include "cgnslib.h"
#include "cgnstypes.h"
#endif // WIN32
#include "ArrayInclude.hpp"
#include "dataStructure.hpp"
#include "geometryModel.hpp"
#include "meshElements.hpp"
#include "zoneMesh.hpp"
namespace OpenHurricane {

    /*!
     * \brief The base class of input and output hdf5 file.
     */
    class cgnsIO {

    public:
        enum accessState : short { OPENED, CLOSED };

        /*!brief Mode used for opening the file*/
        enum CGNSMode {
            READ = CG_MODE_READ,     /**< Read only mode.*/
            WRITE = CG_MODE_WRITE,   /**< Write only mode.*/
            MODIFY = CG_MODE_MODIFY, /**< Reading and /or writing is allowed.*/
            UNKNOWN_MODE
        };

        enum CGNSFileType {
            IS_CGNS_ADF_FILE = CG_FILE_ADF,
            IS_CGNS_ADF2_FILE = CG_FILE_ADF2,
            IS_CGNS_HDF5_FILE = CG_FILE_HDF5,
            IS_CGNS_NONE_FILE = CG_FILE_NONE
        };

        enum class CGNSGridLocation {
            CGNS_NULL = GridLocationNull,
            CGNS_UserDefined = GridLocationUserDefined,
            CGNS_Vertex = Vertex,
            CGNS_CellCenter = CellCenter,
            CGNS_FaceCenter = FaceCenter,
            CGNS_IFaceCenter = IFaceCenter,
            CGNS_JFaceCenter = JFaceCenter,
            CGNS_KFaceCenter = KFaceCenter,
            CGNS_EdgeCenter = EdgeCenter
        };

    protected:
        /**\brief File name.*/
        fileName filename_;

        /**\brief CGNS file index number.*/
        int cgnsFN_;

    public:
        /**\brief CGNS file index number.*/
        hur_nodiscard int cgnsFN() const noexcept { return cgnsFN_; }

    protected:
        /**\brief CGNS file mode.*/
        CGNSMode mode_;

        /*! \brief The access state of the file.*/
        accessState openClosed_;

        /*! \brief Set file opened.*/
        inline void setOpened() noexcept;

        /*! \brief Set file opened.*/
        inline void setClosed() noexcept;

        void setMode(const int flag);

        CGNSGridLocation convertGridLocation(const GridLocation_t &gl) const;
        GridLocation_t convertGridLocation(const CGNSGridLocation &gl) const;

    public:
        /**\brief NUll constructor*/
        inline cgnsIO();

        /**\brief Construct from file name.*/
        inline cgnsIO(const fileName &fN);

        /**\brief Construct from file name.*/
        cgnsIO(const fileName &fN, const int flg);

        /**\brief Destructor.*/
        inline ~cgnsIO() noexcept;

        // FIle operation

        inline void open(const fileName &fN);
        void open();

        inline void open(const unsigned int flg);

        hur_nodiscard float cgVersion();

        hur_nodiscard int cgPrecision();

        /*!\brief Return the file type of cgns.*/
        hur_nodiscard int isCGNS();

        void setCGNSFileType(const int ft);
        hur_nodiscard int getCGNSFileType();

        inline void close();

        /*!\brief Return true if the file opened.*/
        hur_nodiscard inline bool opened() const noexcept;

        /*!\brief Return true if the file closed.*/
        hur_nodiscard inline bool closed() const noexcept;

        /*!\brief Return true if the file is read only.*/
        hur_nodiscard inline bool isRead() const noexcept;

        /*!\brief Return true if the file is write only.*/
        hur_nodiscard inline bool isWrite() const noexcept;

        /*!\brief Return true if the file is modify.*/
        hur_nodiscard inline bool isModify() const noexcept;

        hur_nodiscard DataType_t getDataTypeFromTheProgram() const;

        hur_nodiscard cellShapeType::shapeTypes
        checkElementType(const ElementType_t &type) const noexcept;
        hur_nodiscard ElementType_t
        convertElementType(const cellShapeType::shapeTypes &type) const noexcept;

        hur_nodiscard faceShapeType::shapeTypes
        checkFaceElementType(const ElementType_t &type) const noexcept;
                
        hur_nodiscard ElementType_t
        convertFaceElementType(const faceShapeType::shapeTypes &type) const noexcept;

        hur_nodiscard faceBCType::bcTypes checkBCType(const BCType_t &type) const noexcept;
        hur_nodiscard BCType_t convertBCType(const faceBCType::bcTypes &type) const noexcept;

        /**
         * \brief The number of nodes for the given ElementType.
         */
        hur_nodiscard int numNodes(const ElementType_t &ty) const;

        /*!\brief Return true if the zone is unstructured.*/
        hur_nodiscard bool isUnstructuredZone(const int B, const int Z) const;

        // Read

        /*!\brief Read the number of bases from filename_.*/
        hur_nodiscard int readNBases() const;

        /*!\brief Read the name, cell dimension and physical dimesnion of base B*/
        void readBase(const int B, string &baseName, integer &cellDim, integer &phyDim) const;
        /*!\brief Read the cell dimension of base B.*/
        void readCellDim(const int B, integer &cellDim) const;
        /*!\brief Read the cell dimension of base B.*/
        hur_nodiscard inline int readCellDim(const int B) const;

        /*!\brief Read the number of zones in base B from filename_.*/
        hur_nodiscard int readNZones(const int B) const;
        /*!\brief Read the index dimension for zone Z in base B. (For unstructured zones it will be 1)*/
        hur_nodiscard int readIndexDim(const int B, const int Z) const;
        /*!\brief Read the type of zone.*/
        void readZoneType(const int B, const int Z, ZoneType_t *zoneType) const;
        /*!\brief Read the zone info.*/
        void readZone(const int B, const int Z, string &zoneName, cgsize_t *size) const;
        /*!\brief Read unstructured zone info.
         *   The total number of vertex = NVertex + NBoundVertex
         */
        void readZone(const int B, const int Z, string &zoneName, integer &NVertex,
                      integer &NCell3D, integer &NBoundVertex) const;

        /*!\brief Read simulation type from base B.
         * The valid types are CG_NULL, CG_UserDefined, TimeAccurate, and NonTimeAccurate.
         */
        void readSimulationType(const int B, SimulationType_t &st) const;

        hur_nodiscard int readNCoords(const int B, const int Z) const;

        void readCoordInfo(const int B, const int Z, const int C, DataType_t &dt,
                           string &coordName) const;

        void readCoord(const int B, const int Z, const integer NVertex, const DataType_t dt,
                       const string &coordName, realArray &coord) const;
        hur_nodiscard realArray readCoord(const int B, const int Z, const integer NVertex,
                                          const DataType_t dt, const string &coordName) const;
        void readCoord(const int B, const int Z, const integer NVertex, vectorArray &coord) const;

        hur_nodiscard int readNSections(const int B, const int Z) const;
        void readSection(const int B, const int Z, const int S, string &eleSectName,
                         ElementType_t *type, integer &start, integer &end, int *nbndry,
                         int *parentFlag) const;

        hur_nodiscard integer readElementDataSize(const int B, const int Z, const int S) const;

        void readElement(const int B, const int Z, const int S, cgsize_t *elementData,
                         cgsize_t *parentData) const;

        /*!\brief Read the number of boundary conditions.
         * \param[in] B - The index of the base.
         * \param[in] Z - The index of the zone in the B_th base
         * \return Return the number of the boundary conditions.
         */
        hur_nodiscard int readNBndCon(const int B, const int Z) const;

        /*!\brief Read the id of boundary conditions.
         * \param[in] B - The index of the base.
         * \param[in] Z - The index of the zone in the B_th base
         * \param[in] BC - The index of the boundary condition in the B_th base of Z_th zone
         * \return Return the id of boundary conditions.
         */
        hur_nodiscard double readBndConId(const int B, const int Z, const int BC) const;

        /*!\brief Read the grid location of boundary conditions.
         * \param[in] B - The index of the base.
         * \param[in] Z - The index of the zone in the B_th base
         * \param[in] BC - The index of the boundary condition in the B_th base of Z_th zone
         * \return Return the grid location of boundary conditions.
         */
        hur_nodiscard CGNSGridLocation readBndConLocation(const int B, const int Z,
                                                          const int BC) const;

        void readBndConInfo(const int B, const int Z, const int BC, faceZoneList &fzl) const;

        // Grid Connectivity

        hur_nodiscard int readNConns(const int B, const int Z) const;

        void readConnInfo(const int B, const int Z, const int I, string &connectName) const;

        // Family

        hur_nodiscard int readNFamilies(const int B) const;
        hur_nodiscard int readNodeNFamilies() const;

        // Write

        /*!\brief Write base info and return the base index.*/
        int writeBase(const string &baseName, const integer cellDim, const integer phyDim);

        /*!\brief Write base info and return the base index.*/
        inline int writeBase(const string &baseName);

        /*!\brief Write zone info and return the zone index.*/
        int writeZone(const int B, const string &zoneName, cgsize_t *size, ZoneType_t zoneType);
        /*!\brief Write unstructured zone info and return the zone index.*/
        int writeZone(const int B, const string &zoneName, const integer &NVertex,
                      const integer &NCell3D, const integer &NBoundVertex = 0);

        /*!\brief Write simulation type from base B.
         * The valid types are CG_NULL, CG_UserDefined, TimeAccurate, and NonTimeAccurate.
         */
        void writeSimulationType(const int B, const SimulationType_t st);

        /*!\brief Read the grid location of boundary conditions.
         * \param[in] B - The index of the base.
         * \param[in] Z - The index of the zone in the B_th base
         * \param[in] BC - The index of the boundary condition in the B_th base of Z_th zone
         * \return Return the grid location of boundary conditions.
         */
        void writeBndConLocation(const int B, const int Z, const int BC,
                                 const CGNSGridLocation &gl) const;

        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
         *      Read and write cgns DimensionalUnits_t Nodes                     *
        \* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

        //void readUnits();
        void writeUnits(const int B) const;

        // Parsing

        hur_nodiscard integerListList parsingElementConn(const cgsize_t *hur_restrict ele,
                                                         const integer eleSize,
                                                         const integer eleDataSize,
                                                         const ElementType_t &tpptr) const;

        void parsingCellElementConn(const cgsize_t *hur_restrict ele, const integer start,
                                    const integer end, const integer eleDataSize,
                                    const ElementType_t &tpptr, cellList &cells) const;
    };

} // namespace OpenHurricane

#endif // USES_CGNS

#include "cgnsIO.inl"