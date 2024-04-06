/*!
 * \file HurMPI.cpp
 * \brief Main subroutines for the mpi structures.
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

 #include "HurMPI.hpp"


namespace OpenHurricane {
    hur_nodiscard List<List<int>> HurMPI::getProIdOnNodes() {
#ifdef MPI_PARALLEL
        List<List<int>> pidN(node_);
        if (!multiNodes()) {
            pidN[0].resize(getProcSize());
            for (integer ip = 0; ip < getProcSize(); ++ip) {
                pidN[0][ip] = ip;
            }
        } else {
            int *pidNN = new int[getProcSize()];
            MPI_Allgather(&procNode_, 1, MPI_INT, pidNN, 1, MPI_INT, HurWorldComm_);
            List<int> count(node_);
            count = 0;
            for (integer ip = 0; ip < getProcSize(); ++ip) {
                pidN[pidNN[ip]][count[pidNN[ip]]++] = ip;
            }
            delete[] pidNN;
        }
        return pidN;
#else // Not define MPI_PARALLEL
    List<List<int>> pidN(node_);
    pidN[0].resize(1);
    pidN[0][0] = 0;
    return pidN;
#endif // MPI_PARALLEL
    }
 }
