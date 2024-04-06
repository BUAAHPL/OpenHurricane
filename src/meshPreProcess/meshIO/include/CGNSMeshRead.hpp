/*!
 * \file CGNSMeshRead.hpp
 * \brief Headers of the CGNS mesh reading.
 *        The subroutines and functions are in the <i>CGNSMeshRead.cpp</i> file.
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

#include "originMeshRead.hpp"
#ifdef USES_CGNS

namespace OpenHurricane {
    class cgnsIO;

    class CGNSMeshRead : public originMeshRead {
    private:
        // Private data

        hur_nodiscard int readBase(const cgnsIO &cff) const;

        void checkNumBase(const int nBase) const;
        void checkNumZone(const int nZone) const;

        hur_nodiscard int readNumZone(const cgnsIO &cff, const string &baseName,
                                      const int iBase) const;

        void setZoneSize(const cgnsIO &cff, const int nBase);

        void setPoints(vectorArray &coord);
        void setCells(integerListList &cellConn, const integer start, const integer end);

        void setNGlobalFace();
        void setCellOriginalIndex();

        void formingFaces(const integerListListList &faceZoneEleConn,
                          const integer faceTableCapacity);

        void getFaceZoneEleConn(const cgnsIO &cff, integerListListList &faceZoneEleConn);

        void getGridConnectivity(const cgnsIO &cff, const integer nBase, const integer nZone);

    public:
        // Static data

        /*!\brief Type name.*/
        declareClassNames;

        // Constructors

        /*!\brief Null constructor.*/
        CGNSMeshRead();

        /*!\brief Construct from components.*/
        explicit CGNSMeshRead(const fileName &, const int);

        CGNSMeshRead(const std::string &str, const int);

        /**\brief Destructor.*/
        virtual ~CGNSMeshRead() noexcept {}

    protected:
        virtual void reading(string &gridUnit);
    };

} // namespace OpenHurricane

#endif // USES_CGNS