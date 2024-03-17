/*!
 * \file FluentMeshRead.hpp
 * \brief Headers of the fluent mesh reading.
 *        The subroutines and functions are in the <i>FluentMeshRead.cpp</i> file.
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

#include "originMeshRead.hpp"

namespace OpenHurricane {
    class FluentMeshRead : public originMeshRead {
    private:
        // For fluent mesh

        // Section index
        enum sectionIndexFluent {
            XF_COMMENT = 0,
            XF_HEADER = 1,
            XF_DIMENSION = 2,
            XF_NODE = 10,
            XF_PERIODIC_FACE = 18,
            XF_CELL = 12,
            XF_FACE = 13,
            XF_EDGE = 11,
            XF_FACE_TREE = 59,
            XF_CELL_TREE = 58,
            XF_INTERFACE_PARENTS = 61,
            XF_RP_TV = 39, // Zone section
            XF_RP_TV1 = 45 // or Zone section
        };

        std::map<short int, sectionIndexFluent> sectionIndexFluentMap_;

        void createSectionIndexFluentMap();

        void parsingFluent(std::string &);
        void parsingFluentBinary();

        void skipBrackets(std::ifstream &myin, const int nBrackets) const;
        void skipBrackets(std::ifstream &myin, std::string &str, const int nBrackets) const;
        void skipTo(std::ifstream &myin, const char *c) const;

        void parsingFluentPointBinary(std::ifstream &myin, std::string &line, std::string &index,
                                      integer &countNode);
        void parsingFluentCellBinary(std::ifstream &myin, std::string &line, std::string &index,
                                     integer &totalCells);

        void parsingFluentFaceBinary(std::ifstream &myin, std::string &line, std::string &index,
                                     integer &totalFaces);

        void parsingFluenZoneSectionBinary(std::ifstream &myin, std::string &line);

        void parsingFluentPeriodicFaceBinary(std::ifstream &myin, std::string &line,
                                             std::string &index);

        // Parsing section header for nodes, faces or cells section.
        // (1) The format of nodes section header is as follows
        //      (10 (zone-id first-index last-index bcType ND)(
        //     or
        //      (10 (zone-id first-index last-index bcType ND)
        //      (
        //   When zone-id = 0, first-index will be one, last-index will be
        //   the total number of nodes in hexadecimal, bcType is zero, ND is omitted.
        // (2) The format of faces section header is as follows
        //      (13 (zone-id first-index last-index bc-bcType face-bcType))
        //   When zone-id = 0, first-index will be one, last-index will be
        //   the total number of faces in hexadecimal, bc-bcType and face-bcType are omitted.
        // (3) The format of cells section header is as follows
        //      (12 (zone-id first-index last-index bcType element-bcType))
        //   When zone-id = 0, first-index will be one, last-index will be
        //   the total number of cells in hexadecimal, bcType and element-bcType are omitted.
        // Return 0 if it's the total number section.
        // Return 1 if it's the single section.
        short parsingHeaderStrFluent(const std::string &, integer &zoneId, integer &fi, integer &li,
                                     short &type, short &type2, short &beginListCount,
                                     short &endListCount) const;
        void removeBeginSpaceFluent(std::string &str, size_t &firstSpacePos) const;

        // Low efficiency
        void readFaceConnectionFluent(const faceZone &, const std::string &, const integer,
                                      const short);
        void readFaceConnectionNewFluent(const faceZone &, const std::string &, const integer fI,
                                         const integer lI, const short);
        void readMixedCellShapeFluent(const std::string &, const integer, const integer);
        void readNodesCoordinatesFluent(const std::string &, const integer, const integer,
                                        unsigned short);

        bool isFluentMeshHasBinaryFormat() const;

        void countUnreferencedFaces(integer &n, integerList &unrefFL,
                                    integerList &unrerfFinFZ) const;
        void removeUnreferencedFaces();
        void readCFluent();
        void readFluent();

    public:
        declareClassNames;

        /*!\brief Null constructor.*/
        FluentMeshRead();

        /*!\brief Construct from components.*/
        explicit FluentMeshRead(const fileName &, const int);

        FluentMeshRead(const std::string &str, const int);

        /**\brief Destructor.*/
        virtual ~FluentMeshRead() noexcept {}

        virtual void clear() noexcept;

        inline virtual hur_nodiscard bool checkCoupledWall() const noexcept;

    protected:
        virtual void reading(string &gridUnit);
    };

} // namespace OpenHurricane

#include "FluentMeshRead.inl"