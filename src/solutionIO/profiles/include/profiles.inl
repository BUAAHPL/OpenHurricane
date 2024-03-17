/*!
 * \file profiles.inl
 * \brief In-Line subroutines of the <i>profiles.hpp</i> file.
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

template <class Type>
inline void OpenHurricane::profiles::allGatherVList(const List<Type> &send,
                                                    List<Type> &allL) const {
    integerList nSizeL(HurMPI::getProcSize(), Zero);
    integerList displs(HurMPI::getProcSize(), Zero);
    nSizeL[HurMPI::getProcRank()] = send.size();
    HurMPI::allGatherList(nSizeL, HurMPI::getComm());
    for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
        displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
    }

    integer ssi = 0;
    for (integer i = 0; i < nSizeL.size(); ++i) {
        ssi += nSizeL[i];
    }
    allL.resize(ssi);

    HurMPI::Request request;
    HurMPI::iallGatherv(send.data(), send.size(), feature<Type>::MPIType, allL.data(),
                        nSizeL.data(), displs.data(), feature<Type>::MPIType, HurMPI::getComm(),
                        &request);
    HurMPI::wait(&request, MPI_STATUSES_IGNORE);
}