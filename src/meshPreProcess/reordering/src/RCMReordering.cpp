/*!
 * \file RCMReordering.hpp
 * \brief Headers of the  reversed Cuthill-McKee ordering.
 *        The subroutines and functions are in the <i>RCMReordering.cpp</i> file.
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

#include "RCMReordering.hpp"

 namespace OpenHurricane{
	createClassNameStr(RCMReordering,"RCMReordering");
}

namespace OpenHurricane {
     registerObjFty(reordering,RCMReordering,controller);
}

/*!\brief Reversed Cuthill-Mckee ordering.*/
void OpenHurricane::RCMReordering::RCM() {
    iperm_ = 0;
    // Bubble sort procedure to reorder the adjncy
    // array for each node(from lower degree to high degree)
    for (integer n = 0; n < nvtxs_; n++) {
        for (integer i = xadj_[n]; i <= xadj_[n + 1] - 2; i++) {
            for (integer j = xadj_[n + 1] - 1; j > i; j--) {
                integer ni = adjncy_[j];
                integer ideg = xadj_[ni + 1] - xadj_[ni];
                integer nj = adjncy_[j - 1];
                integer jdeg = xadj_[nj + 1] - xadj_[nj];
                if (ideg < jdeg) {
                    adjncy_[j] = nj;
                    adjncy_[j - 1] = ni;
                }
            }
        }
    } // End bubble sort procedure

    // Check bubble sort procedure
    for (integer n = 0; n < nvtxs_; n++) {
        for (integer i = xadj_[n + 1] - 1; i > xadj_[n]; i--) {
            integer ni = adjncy_[i];
            integer ideg = xadj_[ni + 1] - xadj_[ni];
            for (integer j = xadj_[n]; j < i; j++) {
                integer nj = adjncy_[j];
                integer jdeg = xadj_[nj + 1] - xadj_[nj];
                if (ideg < jdeg) {
                    std::string errMsg;
                    errMsg = "Bubble sort error in cell: ";
                    errMsg += std::to_string(n);
                    errorAbortStr(errMsg);
                }
            }
        }
    } // End checking bubble sort procedure

    // Cuthill-McKee ordering

    integerList marker(nvtxs_, OpenHurricane::Zero);

    //20210510 饶思航 改正 iNode初始取值，应取为-1
    //integer inode = 0;
    integer inode = -1;
    integer ldeg;
    integer root = 0;
    while (inode < nvtxs_ - 1) {
        // Finding root node with lowerest degree between nodes that are not marked
        ldeg = nvtxs_;
        for (integer n = 0; n < nvtxs_; n++) {
            if (marker[n] == 0) {
                integer ideg = xadj_[n + 1] - xadj_[n];
                if (ideg < ldeg) {
                    ldeg = ideg;
                    root = n;
                }
            }
        }
        marker[root] = 1;
        inode++;
        iperm_[inode] = root;

        integer istart = inode;
        integer iend = inode;
        //inode++;

        while (iend < nvtxs_ - 1) {
            integer is = iend + 1;
            integer ie = iend;
            integer signal = 0;

            for (integer i = istart; i <= iend; i++) {
                integer n = iperm_[i];
                for (integer j = xadj_[n]; j <= xadj_[n + 1] - 1; j++) {
                    integer nadj = adjncy_[j];
                    if (marker[nadj] == 0) {
                        marker[nadj] = 1;
                        inode++;
                        iperm_[inode] = nadj;
                        ie++;
                        signal++;
                    }
                }
            }
            istart = is;
            iend = ie;

            if (signal == 0) {
                break;
            }
        }
    }

    // Reverse the result
    reverse(iperm_);
    /*for (integer n = 0; n <nvtxs_ / 2; n++)
    {
            integer i = iperm_[n];
            iperm_[n] = iperm_[nvtxs_ - 1 - n];
            iperm_[nvtxs_ - 1 - n] = i;
    }*/

    for (integer n = 0; n < nvtxs_; n++) {
        integer i = iperm_[n];
        perm_[i] = n;
    }
}
