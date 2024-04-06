/*!
 * \file tecplotWriter.hpp
 * \brief Headers of tecplot writer.
 *        The subroutines and functions are in the <i>tecplotWriter.cpp</i> file.
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
#include "dataStructure.hpp"
#include "fileOsstream.hpp"
#include "globalMesh.hpp"
#include "runtimeMesh.hpp"
#include "tecplotFormat.hpp"

namespace OpenHurricane {
    class tecplotWriter {
        static const int intergerValue_;

        static const char magicNumber_[];

        static const float ZONEMARKER_;

        static const float EOHMARKER_;

        static void checkNumVar(const integer nVar1, const stdStringList &namelist);

    public:
        static void writeASCIIToBinary(const char *str, fileOsstream &fos);

        static void writeInit(fileOsstream &fos);

        static void writeFileHeader(fileOsstream &fos, const std::string &title, const integer nVar,
                                    const stdStringList &varName,
                                    const tecplotFormat::FILETYPE fileType = tecplotFormat::FULL);

        static void writeZoneHeader(fileOsstream &fos, const geometryMesh &mesh,
                                    const real solutionTime, const integer nVar);

        static void writeZoneHeader(fileOsstream &fos, const string &zoneName,
                                    const real solutionTime, const integer nVar, const int ip,
                                    const integer nPts, const integer nCls, const integer nFNgh,
                                    const int ZoneType, const int StrandID = -1);

        static void writeDataSection(fileOsstream &fos, const geometryMesh &mesh,
                                     const integer nVar);

        static void writeDataSectionByMaster(fileOsstream &fos, const geometryMesh &mesh,
                                             const integer nVar);

        static void writeDataSection(fileOsstream &fos, const geometryMesh &mesh,
                                     const integer nVar, const integer fzId);

        template <class Type>
        static void writeData(fileOsstream &fos, const Array<Type> &data, const integer size);

        template <class Type>
        static void writeDataWithFactor(fileOsstream &fos, const Array<Type> &data,
                                        const integer size,
                                        const typename feature<Type>::elementType factor);

        static void writeCell(fileOsstream &fos, const geometryMesh &mesh);

        static void writeCellByMaster(fileOsstream &fos, const geometryMesh &mesh);

        static void writeFaceConnectTecplot(fileOsstream &fos, const geometryMesh &mesh);

        static void writeMesh(fileOsstream &fos, const geometryMesh &mesh);

        static void faceConnectList(integerList &faceMap, const geometryMesh &mesh,
                                    const integer fzId);

        static void writeZoneMesh(fileOsstream &fos, const geometryMesh &mesh, const integer fzId);

        static void writeEOHMARKER(fileOsstream &fos);

        /*!\brief Write results of default variables list to tecplot file.
         *        The call must be collective.
         * \note
         *       Cannot work in multi nodes with MPI.
         * \param[in] fos - file out stream.
         * \param[in] mesh - The geometry mesh.
         */
        static void writeResultToTecplot(fileOsstream &fos, const geometryMesh &mesh);

        /*!\brief Write results of the variables specified by <i>outVarName</i> to tecplot file.
         *        The call must be collective.
         * \note
         *       Cannot work in multi nodes with MPI.
         * \param[in] fos - file out stream.
         * \param[in] mesh - The geometry mesh.
         * \param[in] outVarName - The variables name list for output.
         */
        static void writeResultToTecplot(fileOsstream &fos, const geometryMesh &mesh,
                                         const stringList &outVarName);

        /*!\brief Write results of default variables list to tecplot file by the master.
         *        The call must be collective.
         * \param[in] fos - file out stream.
         * \param[in] mesh - The geometry mesh.
         */
        static void writeResultToTecplotByMaster(fileOsstream &fos, const geometryMesh &mesh);

        /*!\brief Write results of default variables list to tecplot file by the master.
         *        The call must be collective.
         * \param[in] fos - file out stream.
         * \param[in] mesh - The geometry mesh.
         */
        static void writeResultToTecplotByMaster(fileOsstream &fos, const geometryMesh &mesh,
                                                 const stringList &outVarName);

    private:
        static void writeMeshPointsToTecplotByMaster(fileOsstream &fos, const geometryMesh &mesh);
    };

    template <class Type>
    inline void tecplotWriter::writeData(fileOsstream &fos, const Array<Type> &data,
                                         const integer size) {
        integer minSize = min(size, data.size());
        for (int j = 0; j < feature<Type>::nElements_; j++) {
            Array<real> components = data.component(j);
            fos.write(reinterpret_cast<const char *>(&components[0]), minSize * sizeof(real));
        }
    }
    template <class Type>
    inline void
    tecplotWriter::writeDataWithFactor(fileOsstream &fos, const Array<Type> &data,
                                       const integer size,
                                       const typename feature<Type>::elementType factor) {
        integer minSize = min(size, data.size());
        for (int j = 0; j < feature<Type>::nElements_; j++) {
            Array<real> components = data.component(j) * factor;
            fos.write(reinterpret_cast<const char *>(&components[0]), minSize * sizeof(real));
        }
    }
} // namespace OpenHurricane