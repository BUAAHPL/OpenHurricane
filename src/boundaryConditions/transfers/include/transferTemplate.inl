#include "transferTemplate.hpp"
/*!
 * \file transferTemplate.inl
 * \brief In-Line subroutines of the <i>transferTemplate.hpp</i> file.
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

template <class Type>
inline OpenHurricane::processTransfer<Type>::processTransfer(const runtimeMesh &mesh, Array<Type> &var,
                                                         const bool isBlock,
                                                         const bool onlyFirstLayer)
    : mesh_(mesh), var_(var), isBlock_(isBlock), onlyFirstLayer_(onlyFirstLayer), requestCut_(),
      cutPRecvs_(), cutPSends_(), requestPer_(), perPRecvs_(), perPSends_() {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
               "total cells (including ghost cells)");
    }
    const integer ghostLayer = mesh_.ghostCellsType();
    const auto &cutZones_ = mesh_.cutZones();
    const auto &processCut = mesh_.processCut();
    integer czs = cutZones_.size() / ghostLayer;

    if ((czs == 0 || ghostLayer == 0) && cutZones_.size() != 0) {
        LFatal("The size of cutZone is %d, but the layer of ghost cells is %d", cutZones_.size(),
               ghostLayer);
    }
}

template <class Type> inline OpenHurricane::processTransfer<Type>::~processTransfer() noexcept {}

template <class Type> inline void OpenHurricane::processTransfer<Type>::clear() noexcept {
    requestCut_.clear();
    cutPRecvs_.clear();
    cutPSends_.clear();

    requestPer_.clear();
    perPRecvs_.clear();
    perPSends_.clear();
}

template <class Type> inline void OpenHurricane::processTransfer<Type>::transferInit() {
    if (isBlock_) {
        return;
    }
    clear();
    cutNonBlockingTransferInit();
    perNonBlockingTransferInit();
}

template <class Type>
inline void OpenHurricane::processTransfer<Type>::transferInit(const integer layerI) {
    if (isBlock_) {
        return;
    } else {
        clear();
        cutNonBlockingTransferInit(layerI);
        perNonBlockingTransferInit(layerI);
    }
}

template <class Type> inline void OpenHurricane::processTransfer<Type>::transferring() {
    if (isBlock_) {
        transferringBlock();
    } else {
        cutNonBlockingTransferring();
        perNonBlockingTransferring();
    }
}

template <class Type>
inline void OpenHurricane::processTransfer<Type>::transferring(const integer layerI) {
    if (isBlock_) {
        transferringBlock(layerI);
    } else {
        cutNonBlockingTransferring();
        perNonBlockingTransferring();
    }
}

template <class Type> inline void OpenHurricane::processTransfer<Type>::waitAll() {
    HurMPI::Status status;
    HurMPI::waitall(requestCut_.size(), requestCut_.data(), &status);
}

template <class Type> inline void OpenHurricane::processTransfer<Type>::transferringBlock() {
    const integer ghostLayer = mesh_.ghostCellsType();
    const auto &cutZones_ = mesh_.cutZones();
    integer czs = cutZones_.size() / ghostLayer;
    if (HurMPI::parRun()) {
        for (integer layerI = 0; layerI < mesh_.ghostCellsType(); layerI++) {
            if (onlyFirstLayer_ && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                cutTransferringBlock(cutZones_[i]);
            }
        }
    }

    const auto &perZones_ = mesh_.perZones();
    czs = perZones_.size() / ghostLayer;

    for (integer layerI = 0; layerI < mesh_.ghostCellsType(); layerI++) {
        if (onlyFirstLayer_ && layerI >= 1) {
            break;
        }
        for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
            perTransferringBlock(perZones_[i]);
        }
    }
}

template <class Type>
inline void OpenHurricane::processTransfer<Type>::transferringBlock(const integer layerI) {
    const integer ghostLayer = mesh_.ghostCellsType();
    const auto &cutZones_ = mesh_.cutZones();
    integer czs = cutZones_.size() / ghostLayer;

    if (HurMPI::parRun()) {
        if (!onlyFirstLayer_ || layerI < 1) {
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                cutTransferringBlock(cutZones_[i]);
            }
        }
    }

    const auto &perZones_ = mesh_.perZones();
    czs = perZones_.size() / ghostLayer;

    if (!onlyFirstLayer_ || layerI < 1) {
        for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
            perTransferringBlock(perZones_[i]);
        }
    }
}

template <class Type>
inline void OpenHurricane::processTransfer<Type>::cutTransferringBlock(const cutZone &cuts) {
    LFatal("Attempt to access null function");
}

template <class Type>
inline void OpenHurricane::processTransfer<Type>::perTransferringBlock(const perZone &pers) {
    LFatal("Attempt to access null function");
}

template <class Type> inline void OpenHurricane::processTransfer<Type>::cutNonBlockingTransferInit() {
    if (!HurMPI::parRun()) {
        return;
    }

    const auto &cutZones = mesh_.cutZones();
    const auto &processCut = mesh_.processCut();
    const integer ghostLayer = mesh_.ghostCellsType();
    integer czs = cutZones.size() / ghostLayer;
    if ((czs == 0 || ghostLayer == 0) && cutZones.size() != 0) {
        LFatal("The size of cutZone is %d, but the layer of ghost cells is %d", cutZones.size(),
               ghostLayer);
    }
    integer sizeSend = Zero;
    integer sizeRecv = Zero;
    if (onlyFirstLayer_) {
        sizeSend += processCut.zoneIdSend()[0].size();
        sizeRecv += processCut.zoneIdRecv()[0].size();
    } else {
        for (integer layerI = 0; layerI < mesh_.ghostCellsType(); layerI++) {
            sizeSend += processCut.zoneIdSend()[layerI].size();
            sizeRecv += processCut.zoneIdRecv()[layerI].size();
        }
    }

    requestCut_.resize(sizeSend + sizeRecv);

    cutPRecvs_.resize(sizeRecv);
    cutPSends_.resize(sizeSend);

    integer cnt = 0;
    for (integer layerI = 0; layerI < mesh_.ghostCellsType(); layerI++) {
        if (onlyFirstLayer_ && layerI >= 1) {
            break;
        }
        const auto &cutRecv = processCut.zoneIdRecv()[layerI];
        for (integer i = 0; i < cutRecv.size(); ++i) {
            const auto zid = cutRecv[i];
            cutPRecvs_.set(i, new cutProcessReceiv<Type>(cutZones[zid]));
            cutPRecvs_[i].recvNonBlock(&requestCut_[cnt++]);
        }
    }
    for (integer layerI = 0; layerI < mesh_.ghostCellsType(); layerI++) {
        if (onlyFirstLayer_ && layerI >= 1) {
            break;
        }
        const auto &cutSend = processCut.zoneIdSend()[layerI];
        for (integer i = 0; i < cutSend.size(); ++i) {
            const auto zid = cutSend[i];
            cutPSends_.set(i, new cutProcessSend<Type>(cutZones[zid]));
            cutPSends_[i].writeBUf(var_);
            cutPSends_[i].sendNonBlock(&requestCut_[cnt++]);
        }
    }
}

template <class Type>
inline void OpenHurricane::processTransfer<Type>::cutNonBlockingTransferInit(const integer layerI) {
    if (!HurMPI::parRun()) {
        return;
    }
    if (onlyFirstLayer_ && layerI >= 1) {
        return;
    }

    const auto &cutZones = mesh_.cutZones();
    const auto &processCut = mesh_.processCut();
    const integer ghostLayer = mesh_.ghostCellsType();
    integer czs = cutZones.size() / ghostLayer;
    if ((czs == 0 || ghostLayer == 0) && cutZones.size() != 0) {
        LFatal("The size of cutZone is %d, but the layer of ghost cells is %d", cutZones.size(),
               ghostLayer);
    }
    integer sizeSend = Zero;
    integer sizeRecv = Zero;

    sizeSend += processCut.zoneIdSend()[layerI].size();
    sizeRecv += processCut.zoneIdRecv()[layerI].size();

    requestCut_.resize(sizeSend + sizeRecv);

    cutPRecvs_.resize(sizeRecv);
    cutPSends_.resize(sizeSend);

    integer cnt = 0;

    const auto &cutRecv = processCut.zoneIdRecv()[layerI];
    for (integer i = 0; i < cutRecv.size(); ++i) {
        const auto zid = cutRecv[i];
        cutPRecvs_.set(i, new cutProcessReceiv<Type>(cutZones[zid]));
        cutPRecvs_[i].recvNonBlock(&requestCut_[cnt++]);
    }

    const auto &cutSend = processCut.zoneIdSend()[layerI];
    for (integer i = 0; i < cutSend.size(); ++i) {
        const auto zid = cutSend[i];
        cutPSends_.set(i, new cutProcessSend<Type>(cutZones[zid]));
        cutPSends_[i].writeBUf(var_);
        cutPSends_[i].sendNonBlock(&requestCut_[cnt++]);
    }
}

template <class Type> inline void OpenHurricane::processTransfer<Type>::perNonBlockingTransferInit() {
    const integer ghostLayer = mesh_.ghostCellsType();
    const auto &perZones = mesh_.perZones();
    if (perZones.size() == 0) {
        return;
    }
    integer czs = perZones.size() / ghostLayer;
    if ((czs == 0 || ghostLayer == 0) && perZones.size() != 0) {
        LFatal("The size of perZone is %d, but the layer of ghost cells is %d", perZones.size(),
               ghostLayer);
    }

    for (integer layerI = 0; layerI < mesh_.ghostCellsType(); layerI++) {
        if (onlyFirstLayer_ && layerI >= 1) {
            break;
        }
        for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
            if (perZones[i].sendProc() == perZones[i].receivProc()) {
                perTransferringBlock(perZones[i]);
            }
        }
    }

    if (!HurMPI::parRun()) {
        return;
    }

    const auto &processPer = mesh_.processPer();
    integer sizeSend = Zero;
    integer sizeRecv = Zero;
    if (onlyFirstLayer_) {
        sizeSend += processPer.zoneIdSend()[0].size();
        sizeRecv += processPer.zoneIdRecv()[0].size();
    } else {
        for (integer layerI = 0; layerI < mesh_.ghostCellsType(); layerI++) {
            sizeSend += processPer.zoneIdSend()[layerI].size();
            sizeRecv += processPer.zoneIdRecv()[layerI].size();
        }
    }

    requestPer_.resize(sizeSend + sizeRecv);

    perPRecvs_.resize(sizeRecv);
    perPSends_.resize(sizeSend);

    integer cnt = 0;
    for (integer layerI = 0; layerI < mesh_.ghostCellsType(); layerI++) {
        if (onlyFirstLayer_ && layerI >= 1) {
            break;
        }
        const auto &perRecv = processPer.zoneIdRecv()[layerI];
        for (integer i = 0; i < perRecv.size(); ++i) {
            const auto zid = perRecv[i];
            perPRecvs_.set(i, new perProcessReceiv<Type>(perZones[zid]));
            perPRecvs_[i].recvNonBlock(&requestPer_[cnt++]);
        }
    }
    for (integer layerI = 0; layerI < mesh_.ghostCellsType(); layerI++) {
        if (onlyFirstLayer_ && layerI >= 1) {
            break;
        }
        const auto &perSend = processPer.zoneIdSend()[layerI];
        for (integer i = 0; i < perSend.size(); ++i) {
            const auto zid = perSend[i];
            perPSends_.set(i, new perProcessSend<Type>(perZones[zid]));
            perPSends_[i].writeBUf(var_);
            perPSends_[i].sendNonBlock(&requestPer_[cnt++]);
        }
    }
}

template <class Type>
inline void OpenHurricane::processTransfer<Type>::perNonBlockingTransferInit(const integer layerI) {
    const integer ghostLayer = mesh_.ghostCellsType();
    const auto &perZones = mesh_.perZones();
    integer czs = perZones.size() / ghostLayer;

    if (onlyFirstLayer_ && layerI >= 1) {
        return;
    }
    for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
        if (perZones[i].sendProc() == perZones[i].receivProc()) {
            perTransferringBlock(perZones[i]);
        }
    }

    if (!HurMPI::parRun()) {
        return;
    }

    const auto &processPer = mesh_.processPer();
    integer sizeSend = Zero;
    integer sizeRecv = Zero;

    sizeSend += processPer.zoneIdSend()[layerI].size();
    sizeRecv += processPer.zoneIdRecv()[layerI].size();

    requestPer_.resize(sizeSend + sizeRecv);

    perPRecvs_.resize(sizeRecv);
    perPSends_.resize(sizeSend);

    integer cnt = 0;

    const auto &perRecv = processPer.zoneIdRecv()[layerI];
    for (integer i = 0; i < perRecv.size(); ++i) {
        const auto zid = perRecv[i];
        perPRecvs_.set(i, new perProcessReceiv<Type>(perZones[zid]));
        perPRecvs_[i].recvNonBlock(&requestPer_[cnt++]);
    }

    const auto &perSend = processPer.zoneIdSend()[layerI];
    for (integer i = 0; i < perSend.size(); ++i) {
        const auto zid = perSend[i];
        perPSends_.set(i, new perProcessSend<Type>(perZones[zid]));
        perPSends_[i].writeBUf(var_);
        perPSends_[i].sendNonBlock(&requestPer_[cnt++]);
    }
}

template <class Type> inline void OpenHurricane::processTransfer<Type>::cutNonBlockingTransferring() {
    if (!HurMPI::parRun() || requestCut_.size() == 0) {
        return;
    }

    HurMPI::waitall(requestCut_.size(), requestCut_.data(), MPI_STATUSES_IGNORE);

    for (integer i = 0; i < cutPRecvs_.size(); ++i) {
        cutPRecvs_[i].readBUf(var_);
    }
}

template <class Type> inline void OpenHurricane::processTransfer<Type>::perNonBlockingTransferring() {
    if (!HurMPI::parRun() || requestPer_.size() == 0) {
        return;
    }

    HurMPI::waitall(requestPer_.size(), requestPer_.data(), MPI_STATUSES_IGNORE);

    for (integer i = 0; i < perPRecvs_.size(); ++i) {
        perPRecvs_[i].readBUf(var_);
    }
}
