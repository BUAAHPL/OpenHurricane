/*!
 * \file HurricaneMeshRead.hpp
 * \brief Headers of the Hurricane mesh reading.
 *        The subroutines and functions are in the <i>HurricaneMeshRead.cpp</i> file.
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
    class HurricaneMeshRead : public originMeshRead {
    private:
        // Private functions

        void readOriginAndAxis(const hdf5I &fos, const string &gridGroupName);
        void readCellZonesHurricane(const hdf5I &fos, const string &gridGroupName,
                                    const integer nCZ);
        void setCellsHurricane(const hdf5I &fos, const string &gridGroupName, const integer nC);
        void readFaceZonesHurricane(const hdf5I &fos, const string &gridGroupName,
                                    const integer nFZ);
        void setFacesHurricane(const hdf5I &fos, const string &gridGroupName, const integer nF);
        void readFacesHurricane(const hdf5I &fos, const string &gridGroupName, const faceZone &fz,
                                const integerArrayArray &faceConnect);
        void readPointZonesHurricane(const hdf5I &fos, const string &gridGroupName,
                                     const integer nPZ);
        void setPointsHurricane(const hdf5I &fos, const string &gridGroupName, const integer nP);
        void readPointsHurricane(const hdf5I &fos, const string &gridGroupName, const pointZone &pz,
                                 const vectorArray &tp);
        void setPeriodicPairListSizeHurricane(const hdf5I &fos, const string &gridGroupName,
                                              const integer nPP);
        void readPeriodicPairListHurricane(const hdf5I &fos, const string &gridGroupName);

        integer readGridIntAttrHurricane(const hdf5I &fos, const string &gridGroupName,
                                         const string &attrName);

        void readHurricane(const hdf5I &fos, const string &gridGroupName, string &gridUnit);

    public:
        declareClassNames;

        /*!\brief Null constructor.*/
        HurricaneMeshRead();

        /*!\brief Construct from components.*/
        explicit HurricaneMeshRead(const fileName &, const int);

        HurricaneMeshRead(const std::string &str, const int);

        /**\brief Destructor.*/
        virtual ~HurricaneMeshRead() noexcept {}

        inline virtual hur_nodiscard bool isHurricaneMesh() const noexcept;

    protected:
        virtual void reading(string &gridUnit);

    protected:

        /**
         * \brief Read the loas weight from mesh file.
         * \param[in] fos - File operator.
         * \param[in] gridGroupName - The grid group name in file fos.
         * \param[in] nC - The number of cells.
         * \param[in] nCZ - The number of cell zones.
         */
        void readCellLoadWeightHurricane(const hdf5I &fos, const string &gridGroupName,
                                         const integer nC, const integer nCZ);
    };

} // namespace OpenHurricane

#include "HurricaneMeshRead.inl"