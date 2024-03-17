/*!
 * \file CRSMatrixProcessTransfer.cpp
 * \brief Main subroutines of CRSMatrixProcessTransfer.
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

#include "CRSMatrixProcessTransfer.hpp"
template <>
void OpenHurricane::CRSMatrixProcessTransfer<OpenHurricane::realArray>::transferData(
    const Array<realArray> &var, Array<realArray> &recvVar) {
    if (!HurMPIBase::parRun()) {
        for (integer i = 0; i < perSameProc_.size(); ++i) {
            perSameProc_[i].transferData(var, recvVar);
        }
    } else {

        integer cntCutRe = 0;
        for (integer i = 0; i < cutPRecvs_.size(); ++i) {
            cutPRecvs_[i].recvNonBlock(&requestCut_[cntCutRe++],
                                       cutPRecvs_[i].zone().des().size() * var.first().size());
        }

        for (integer i = 0; i < cutPSends_.size(); ++i) {
            cutPSends_[i].writeBUf(var);
            cutPSends_[i].sendNonBlock(&requestCut_[cntCutRe++]);
        }

        integer cntPerRe = 0;
        for (integer i = 0; i < perPRecvs_.size(); ++i) {
            perPRecvs_[i].recvNonBlock(&requestPer_[cntPerRe++],
                                       perPRecvs_[i].zone().des().size() * var.first().size());
        }

        for (integer i = 0; i < perPSends_.size(); ++i) {
            perPSends_[i].writeBUf(var);
            perPSends_[i].sendNonBlock(&requestPer_[cntPerRe++]);
        }

        for (integer i = 0; i < perSameProc_.size(); ++i) {
            perSameProc_[i].transferData(var, recvVar);
        }

        HurMPIBase::waitall(requestCut_.size(), requestCut_.data(), MPI_STATUSES_IGNORE);

        for (integer i = 0; i < cutPRecvs_.size(); ++i) {
            cutPRecvs_[i].readBUf(recvVar);
        }

        HurMPIBase::waitall(requestPer_.size(), requestPer_.data(), MPI_STATUSES_IGNORE);
        for (integer i = 0; i < perPRecvs_.size(); ++i) {
            perPRecvs_[i].readBUf(recvVar);
        }
    }
}