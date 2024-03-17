/*!
 * \file multilevelNestedDissectionByMetis.hpp
 * \brief Headers of the multilevel nested dissection algorithm by metis.
 *        The subroutines and functions are in the <i>multilevelNestedDissectionByMetis.cpp</i> file.
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

#include "multilevelNestedDissectionByMetis.hpp"
#include "logFile.hpp"
#include "metis.h"

 namespace OpenHurricane{
	createClassNameStr(multilevelNestedDissectionByMetis,"MetisReordering");
}

namespace OpenHurricane {
     registerObjFty(reordering,multilevelNestedDissectionByMetis,controller);
}

void OpenHurricane::multilevelNestedDissectionByMetis::reorder() {
    Pout << "    Info: reordering mesh using the multilevel nested dissection "
            "algorithm (METIS_NodeND)..."
         << std::endl;
    integer nb1 = getBandwidth();
    computeArray();

    idx_t options[METIS_NOPTIONS];

    METIS_SetDefaultOptions(options);

    // options[METIS OPTION NUMBERING]
    //	 Used to indicate which numbering scheme is used for the adjacency structure of a graph or the elementnode
    //	 structure of a mesh.The possible values are :
    //    0 C - style numbering is assumed that starts from 0.
    //	  1 Fortran - style numbering is assumed that starts from 1.
    options[METIS_OPTION_NUMBERING] = 0;

    options[METIS_OPTION_NITER] = 10;

    idx_t *vwgt;

    vwgt = nullptr;

    int returnFlag;

    returnFlag = METIS_NodeND(&nvtxs_, xadj_.data(), adjncy_.data(), vwgt, options, perm_.data(),
                              iperm_.data());

    if (returnFlag == METIS_ERROR_INPUT) {
        LFatal("An input error occur in metis function!");

    } else if (returnFlag == METIS_ERROR_MEMORY) {
        LFatal("Could not allocate the required memory in metis function!");

    } else if (returnFlag == METIS_ERROR) {
        LFatal("Unknown error in metis functions!");
    }

    update();

    integer nb2 = getBandwidth();

    printBandwidthReduction(nb1, nb2);
}
