﻿/*!
 * \file relayScatter.cpp
 * \brief Main subroutines for relayScatter.
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
#include "relayScatter.hpp"

template <>
void OpenHurricane::relayScatterFunc::relayScatter(const integerArrayArray &cellOfProc,
                                               const Array<real> &allArray, Array<real> &ipArray) {
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
                Array<real> tArray(cellSize[ip]);
                for (integer i = 0; i < cellOfProc[ip].size(); ++i) {
                    tArray[i] = allArray[cellOfProc[ip][i]];
                }
                HurMPI::sendList(tArray, ip, HurMPI::masterNo(), HurMPI::getComm());
            }
        }
    } else {
        HurMPI::Status st;
        HurMPI::recvList(ipArray, HurMPI::masterNo(), HurMPI::masterNo(), HurMPI::getComm(), &st);
    }
}
