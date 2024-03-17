/*!
 * \file relayScatter.inl
 * \brief In-Line subroutines of the <i>relayScatter.hpp</i> file.
 * \author Yang Hongzhen
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
#include "relayScatter.hpp"

template <class Type>
inline void OpenHurricane::relayScatterFunc::relayScatter(const integerArrayArray &cellOfProc,
                                                      const Array<Type> &allArray,
                                                      Array<Type> &ipArray) {
    integerList cellSize(HurMPI::getProcSize(), Zero);
    if (HurMPI::master()) {
        for (integer ip = 0; ip < cellOfProc.size(); ++ip) {
            cellSize[ip] = cellOfProc[ip].size();
        }
    }
    HurMPI::bcastList(cellSize, HurMPI::masterNo(), HurMPI::getComm());
    ipArray.resize(cellSize[HurMPI::getProcRank()]);

    if (HurMPI::master()) {
        for (integer ip = 0; ip < cellOfProc.size(); ++ip) {
            if (HurMPI::masterNo() == ip) {
                for (integer i = 0; i < cellOfProc[ip].size(); ++i) {
                    ipArray[i] = allArray[cellOfProc[ip][i]];
                }
            } else {
                for (integer k = 0; k < feature<Type>::nElements_; ++k) {
                    Array<typename feature<Type>::elementType> tArray(cellSize[ip]);
                    for (integer i = 0; i < cellOfProc[ip].size(); ++i) {
                        tArray[i] = allArray[cellOfProc[ip][i]][k];
                    }
                    HurMPI::sendList(tArray, ip, k, HurMPI::getComm());
                }
            }
        }
    } else {
        for (integer k = 0; k < feature<Type>::nElements_; ++k) {
            Array<typename feature<Type>::elementType> tArray(cellSize[HurMPI::getProcRank()]);
            HurMPI::Status st;
            HurMPI::recvList(tArray, HurMPI::masterNo(), k, HurMPI::getComm(), &st);
            for (integer i = 0; i < tArray.size(); ++i) {
                ipArray[i][k] = tArray[i];
            }
        }
    }
}

template <class Type>
inline void
OpenHurricane::relayScatterFunc::relayReorder(const Array<Type> &pArray, Array<Type> &ipArray,
                                          const integerList &perm_, const integerList &iperm_) {
    for (integer n = 0; n < pArray.size(); n++) {
        integer i = iperm_[n];
        ipArray[n] = pArray[i];
    }
}