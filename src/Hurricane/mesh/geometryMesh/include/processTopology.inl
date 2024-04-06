#include "processTopology.hpp"
/*!
 * \file processTopology.inl
 * \brief In-Line subroutines of the <i>processTopology.hpp</i> file.
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

template <class zoneList>
inline void OpenHurricane::processTopology<zoneList>::getZoneId(const integer nLayer) {
    if (!HurMPI::parRun()) {
        return;
    }
    zoneIdSend_.resize(nLayer);
    zoneIdRecv_.resize(nLayer);

    integer czs = zones_.size() / nLayer;
    if (czs == 0 && zones_.size() != 0) {
        LFatal("The size of zone is %d, but the layer of ghost cells is %d", zones_.size(), nLayer);
    }
    integerList countSend(nLayer, Zero);
    integerList countRecv(nLayer, Zero);
    for (integer layerI = 0; layerI < nLayer; layerI++) {
        for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
            if (HurMPI::isThisProc(zones_[i].sendProc()) &&
                !HurMPI::isThisProc(zones_[i].receivProc())) {
                countSend[layerI]++;
            }
            if (!HurMPI::isThisProc(zones_[i].sendProc()) &&
                HurMPI::isThisProc(zones_[i].receivProc())) {
                countRecv[layerI]++;
            }
        }
    }

    for (integer layerI = 0; layerI < nLayer; layerI++) {
        zoneIdSend_[layerI].resize(countSend[layerI], -1);
        zoneIdRecv_[layerI].resize(countRecv[layerI], -1);

        countSend[layerI] = 0;
        countRecv[layerI] = 0;
    }

    for (integer layerI = 0; layerI < nLayer; layerI++) {
        for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
            if (HurMPI::isThisProc(zones_[i].sendProc()) &&
                !HurMPI::isThisProc(zones_[i].receivProc())) {
                zoneIdSend_[layerI][countSend[layerI]] = i;
                countSend[layerI]++;
            }
            if (!HurMPI::isThisProc(zones_[i].sendProc()) &&
                HurMPI::isThisProc(zones_[i].receivProc())) {
                zoneIdRecv_[layerI][countRecv[layerI]] = i;
                countRecv[layerI]++;
            }
        }
    }
}

template <class zoneList>
inline OpenHurricane::processTopology<zoneList>::processTopology(const zoneList &zones,
                                                             const integer nlayer)
    : zones_(zones), zoneIdSend_(), zoneIdRecv_() {
    getZoneId(nlayer);
}

template <class zoneList>
inline OpenHurricane::processTopology<zoneList>::~processTopology() noexcept {}

template <class zoneList>
hur_nodiscard inline const zoneList &OpenHurricane::processTopology<zoneList>::zones() const noexcept {
    return zones_;
}

template <class zoneList>
hur_nodiscard inline const OpenHurricane::integerListList &
OpenHurricane::processTopology<zoneList>::zoneIdSend() const noexcept {
    return zoneIdSend_;
}

template <class zoneList>
hur_nodiscard inline const OpenHurricane::integerListList &
OpenHurricane::processTopology<zoneList>::zoneIdRecv() const noexcept {
    return zoneIdRecv_;
}

template <class zoneList>
hur_nodiscard inline OpenHurricane::integer
OpenHurricane::processTopology<zoneList>::zoneIdSend(const integer nLayer, const integer n) const {
    if (!HurMPI::parRun()) {
        return -1;
    }
    return zoneIdSend_[nLayer][n];
}

template <class zoneList>
hur_nodiscard inline OpenHurricane::integer
OpenHurricane::processTopology<zoneList>::zoneIdRecv(const integer nLayer, const integer n) const {
    if (!HurMPI::parRun()) {
        return -1;
    }
    return zoneIdRecv_[nLayer][n];
}

template <class zoneType, class dataType>
inline OpenHurricane::processSend<zoneType, dataType>::processSend(const zoneType &zone)
    : zone_(zone) {}

template <class zoneType, class dataType>
inline OpenHurricane::processSend<zoneType, dataType>::~processSend() noexcept {}

template <class zoneType, class dataType>
inline void OpenHurricane::processSend<zoneType, dataType>::clear() noexcept {
    buf_.clear();
}

template <class zoneType, class dataType>
inline void OpenHurricane::processSend<zoneType, dataType>::writeBUf(const Array<dataType> &f) {
    buf_.resize(zone_.sor().size() * feature<dataType>::nElements_, Zero);

    for (integer i = 0; i < zone_.sor().size(); i++) {
        for (int j = 0; j < feature<dataType>::nElements_; j++) {
            buf_[i * feature<dataType>::nElements_ + j] = f[zone_.sor()[i]][j];
        }
    }
}

template <class zoneType, class dataType>
inline void OpenHurricane::processSend<zoneType, dataType>::send() const {
    HurMPI::send(buf_.data(), buf_.size(),
                 feature<typename feature<dataType>::elementType>::MPIType, zone_.receivProc(),
                 zone_.sendProc(), HurMPI::getComm());
}

template <class zoneType, class dataType>
inline void OpenHurricane::processSend<zoneType, dataType>::sendNonBlock(HurMPI::Request *request) {
    auto ierr = HurMPI::isend(buf_.data(), buf_.size(),
                              feature<typename feature<dataType>::elementType>::MPIType,
                              zone_.receivProc(), zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Send data error: %s", estr.c_str());
    }
}

template <class zoneType, class dataType>
inline OpenHurricane::processReceiv<zoneType, dataType>::processReceiv(const zoneType &zone)
    : zone_(zone) {}

template <class zoneType, class dataType>
inline OpenHurricane::processReceiv<zoneType, dataType>::~processReceiv() noexcept {}

template <class zoneType, class dataType>
inline void OpenHurricane::processReceiv<zoneType, dataType>::clear() noexcept {
    buf_.clear();
}

template <class zoneType, class dataType>
inline void OpenHurricane::processReceiv<zoneType, dataType>::readBUf(Array<dataType> &f) {
    if (buf_.size() < zone_.des().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        for (int j = 0; j < feature<dataType>::nElements_; j++) {
            f[zone_.des()[i]][j] = buf_[i * feature<dataType>::nElements_ + j];
        }
    }
}

template <class zoneType, class dataType>
inline void OpenHurricane::processReceiv<zoneType, dataType>::recv() {
    buf_.resize(zone_.des().size() * feature<dataType>::nElements_, Zero);
    HurMPI::recv(buf_.data(), buf_.size(),
                 feature<typename feature<dataType>::elementType>::MPIType, zone_.sendProc(),
                 zone_.sendProc(), HurMPI::getComm());
}

template <class zoneType, class dataType>
inline void OpenHurricane::processReceiv<zoneType, dataType>::recvNonBlock(HurMPI::Request *request) {
    buf_.resize(zone_.des().size() * feature<dataType>::nElements_, Zero);
    auto ierr = HurMPI::irecv(buf_.data(), buf_.size(),
                              feature<typename feature<dataType>::elementType>::MPIType,
                              zone_.sendProc(), zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Receiv data error: %s", estr.c_str());
    }
}

template <class zoneType, class dataType>
inline void OpenHurricane::processReceiv<zoneType, dataType>::recvNonBlock(HurMPI::Request *request,
                                                                       const integer sizeBuf) {
    buf_.resize(sizeBuf, Zero);
    auto ierr = HurMPI::irecv(buf_.data(), buf_.size(),
                              feature<typename feature<dataType>::elementType>::MPIType,
                              zone_.sendProc(), zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Receiv data error: %s", estr.c_str());
    }
}

template <class zoneType, class dataType>
hur_nodiscard inline OpenHurricane::Array<typename OpenHurricane::feature<dataType>::elementType> &
OpenHurricane::processReceiv<zoneType, dataType>::buf() noexcept {
    return buf_;
}