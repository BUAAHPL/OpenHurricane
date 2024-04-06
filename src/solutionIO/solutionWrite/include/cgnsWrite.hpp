/*!
 * \file cgnsWrite.hpp
 * \brief Headers of class of cgns file write-out.
 *        The subroutines and functions are in the <i>cgnsWrite.cpp</i> file.
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

#ifdef USES_CGNS
#include "cgnsIO.hpp"
#include "solutionWrite.hpp"

namespace OpenHurricane {
    /*!\brief The base class of cgnsWrite.*/
    class cgnsWrite : public solutionWrite {
    private:
        cgnsIO::CGNSFileType fileType_;

    public:
        declareClassNames;

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        cgnsWrite(const flowModel &flows, const iteration &iter, const runtimeMesh &mesh,
                  const controller &cont);

        /**
         * \brief Destructor.
         */
        virtual ~cgnsWrite() noexcept {}

        virtual void writeToFile() const;

    private:
        mutable bool firstCallCGBC_;
        mutable List<List<integer>> BCParentElem_;
        mutable List<List<integer>> BCParentElemPos_;
        void writePoints(cgnsIO &cff, int B, int Z) const;
        void writeCells(cgnsIO &cff, int B, int Z) const;
        void writeBCFaces(cgnsIO &cff, int B, int Z) const;
        void writeSolutions(cgnsIO &cff, int B, int Z) const;

        int writeArray(cgnsIO &cff, int B, int Z, int S, realArray &tmpv, const char *nstr) const;
    };
} // namespace OpenHurricane

#endif // USES_CGNS