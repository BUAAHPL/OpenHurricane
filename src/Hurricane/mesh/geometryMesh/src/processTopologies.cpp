/*!
 * \file processTopologies.cpp
 * \brief Main subroutines for process topologies.
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

#include "processTopologies.hpp"

//==================================================================================
// For real

template <> void OpenHurricane::cutProcessSend<OpenHurricane::real>::writeBUf(const Array<real> &f) {
    buf_.resize(zone_.sor().size(), Zero);

    for (integer i = 0; i < zone_.sor().size(); i++) {
        buf_[i] = f[zone_.sor()[i]];
    }
}

template <> void OpenHurricane::perProcessSend<OpenHurricane::real>::writeBUf(const Array<real> &f) {
    buf_.resize(zone_.sor().size(), Zero);

    for (integer i = 0; i < zone_.sor().size(); i++) {
        buf_[i] = f[zone_.sor()[i]];
    }
}

//==================================================================================
// For integer

template <> void OpenHurricane::cutProcessSend<OpenHurricane::integer>::writeBUf(const Array<integer> &f) {
    buf_.resize(zone_.sor().size(), Zero);

    for (integer i = 0; i < zone_.sor().size(); i++) {
        buf_[i] = f[zone_.sor()[i]];
    }
}

template <> void OpenHurricane::perProcessSend<OpenHurricane::integer>::writeBUf(const Array<integer> &f) {
    buf_.resize(zone_.sor().size(), Zero);

    for (integer i = 0; i < zone_.sor().size(); i++) {
        buf_[i] = f[zone_.sor()[i]];
    }
}

//==================================================================================
// For complex

template <> void OpenHurricane::cutProcessSend<OpenHurricane::complex>::writeBUf(const Array<complex> &f) {
    buf_.resize(2 * zone_.sor().size(), Zero);

    for (integer i = 0; i < zone_.sor().size(); i++) {
        buf_[2 * i] = f[zone_.sor()[i]].real();
        buf_[2 * i + 1] = f[zone_.sor()[i]].imag();
    }
}

template <> void OpenHurricane::perProcessSend<OpenHurricane::complex>::writeBUf(const Array<complex> &f) {
    buf_.resize(2 * zone_.sor().size(), Zero);

    for (integer i = 0; i < zone_.sor().size(); i++) {
        buf_[2 * i] = f[zone_.sor()[i]].real();
        buf_[2 * i + 1] = f[zone_.sor()[i]].imag();
    }
}

//==================================================================================
// For realArray

template <>
void OpenHurricane::cutProcessSend<OpenHurricane::realArray>::writeBUf(const Array<realArray> &f) {
#ifdef HUR_DEBUG
    auto s0 = f.first().size();
    for (const auto &e : f) {
        if (e.size() != s0) {
            LFatal("Not all element of realArrayArray has the same size.");
        }
    }
#endif // HUR_DEBUG

    buf_.resize(zone_.sor().size() * f.first().size(), Zero);
    for (integer i = 0; i < zone_.sor().size(); i++) {
        for (integer j = 0; j < f.first().size(); ++j) {
            buf_[i * f.first().size() + j] = f[zone_.sor()[i]][j];
        }
    }
}

template <>
void OpenHurricane::perProcessSend<OpenHurricane::realArray>::writeBUf(const Array<realArray> &f) {
#ifdef HUR_DEBUG
    auto s0 = f.first().size();
    for (const auto &e : f) {
        if (e.size() != s0) {
            LFatal("Not all element of realArrayArray has the same size.");
        }
    }
#endif // HUR_DEBUG

    buf_.resize(zone_.sor().size() * f.first().size(), Zero);
    for (integer i = 0; i < zone_.sor().size(); i++) {
        for (integer j = 0; j < f.first().size(); ++j) {
            buf_[i * f.first().size() + j] = f[zone_.sor()[i]][j];
        }
    }
}

//==================================================================================
// For real

template <> void OpenHurricane::cutProcessReceiv<OpenHurricane::real>::readBUf(Array<real> &f) {
    if (buf_.size() != zone_.des().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        f[zone_.des()[i]] = buf_[i];
    }
}

template <> void OpenHurricane::perProcessReceiv<OpenHurricane::real>::readBUf(Array<real> &f) {
    if (buf_.size() != zone_.des().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        f[zone_.des()[i]] = buf_[i];
    }
}

template <>
void OpenHurricane::cutProcessReceiv<OpenHurricane::real>::recvNonBlock(HurMPI::Request *request) {
    buf_.resize(zone_.des().size(), Zero);
    auto ierr = HurMPI::irecv(buf_.data(), buf_.size(), feature<real>::MPIType, zone_.sendProc(),
                              zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Receiv data error: %s", estr.c_str());
    }
}

template <>
void OpenHurricane::perProcessReceiv<OpenHurricane::real>::recvNonBlock(HurMPI::Request *request) {
    if (zone_.sendProc() == zone_.receivProc()) {
        return;
    }

    buf_.resize(zone_.des().size(), Zero);
    auto ierr = HurMPI::irecv(buf_.data(), buf_.size(), feature<real>::MPIType, zone_.sendProc(),
                              zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Receiv data error: %s", estr.c_str());
    }
}

//==================================================================================
// For integer

template <> void OpenHurricane::cutProcessReceiv<OpenHurricane::integer>::readBUf(Array<integer> &f) {
    if (buf_.size() != zone_.des().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        f[zone_.des()[i]] = buf_[i];
    }
}

template <> void OpenHurricane::perProcessReceiv<OpenHurricane::integer>::readBUf(Array<integer> &f) {
    if (buf_.size() != zone_.des().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        f[zone_.des()[i]] = buf_[i];
    }
}

template <>
void OpenHurricane::cutProcessReceiv<OpenHurricane::integer>::recvNonBlock(HurMPI::Request *request) {
    buf_.resize(zone_.des().size(), Zero);
    auto ierr = HurMPI::irecv(buf_.data(), buf_.size(), feature<integer>::MPIType, zone_.sendProc(),
                              zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Receiv data error: %s", estr.c_str());
    }
}

template <>
void OpenHurricane::perProcessReceiv<OpenHurricane::integer>::recvNonBlock(HurMPI::Request *request) {
    if (zone_.sendProc() == zone_.receivProc()) {
        return;
    }

    buf_.resize(zone_.des().size(), Zero);
    auto ierr = HurMPI::irecv(buf_.data(), buf_.size(), feature<integer>::MPIType, zone_.sendProc(),
                              zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Receiv data error: %s", estr.c_str());
    }
}

//==================================================================================
// For complex

template <> void OpenHurricane::cutProcessReceiv<OpenHurricane::complex>::readBUf(Array<complex> &f) {
    if (buf_.size() != 2 * zone_.des().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        f[zone_.des()[i]].real(buf_[2 * i]);
        f[zone_.des()[i]].imag(buf_[2 * i + 1]);
    }
}

template <> void OpenHurricane::perProcessReceiv<OpenHurricane::complex>::readBUf(Array<complex> &f) {
    if (buf_.size() != 2 * zone_.des().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        f[zone_.des()[i]].real(buf_[2 * i]);
        f[zone_.des()[i]].imag(buf_[2 * i + 1]);
    }
}

template <>
void OpenHurricane::cutProcessReceiv<OpenHurricane::complex>::recvNonBlock(HurMPI::Request *request) {
    buf_.resize(2 * zone_.des().size(), Zero);
    auto ierr =
        HurMPI::irecv(buf_.data(), buf_.size(), feature<feature<complex>::elementType>::MPIType,
                      zone_.sendProc(), zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Receiv data error: %s", estr.c_str());
    }
}

template <>
void OpenHurricane::perProcessReceiv<OpenHurricane::complex>::recvNonBlock(HurMPI::Request *request) {
    if (zone_.sendProc() == zone_.receivProc()) {
        return;
    }

    buf_.resize(2 * zone_.des().size(), Zero);
    auto ierr =
        HurMPI::irecv(buf_.data(), buf_.size(), feature<feature<complex>::elementType>::MPIType,
                      zone_.sendProc(), zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Receiv data error: %s", estr.c_str());
    }
}

//==================================================================================
// For realArray

template <> void OpenHurricane::cutProcessReceiv<OpenHurricane::realArray>::readBUf(Array<realArray> &f) {
    if (buf_.size() != zone_.des().size() * f.first().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        for (integer j = 0; j < f.first().size(); ++j) {
            f[zone_.des()[i]][j] = buf_[i * f.first().size() + j];
        }
    }
}

template <> void OpenHurricane::perProcessReceiv<OpenHurricane::realArray>::readBUf(Array<realArray> &f) {
    if (buf_.size() != zone_.des().size() * f.first().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        for (integer j = 0; j < f.first().size(); ++j) {
            f[zone_.des()[i]][j] = buf_[i * f.first().size() + j];
        }
    }
}

template <>
void OpenHurricane::cutProcessReceiv<OpenHurricane::realArray>::recvNonBlock(HurMPI::Request *request) {
    LFatal("This function can not be used for Array<realArray>.");
}

template <>
void OpenHurricane::perProcessReceiv<OpenHurricane::realArray>::recvNonBlock(HurMPI::Request *request) {
    if (zone_.sendProc() == zone_.receivProc()) {
        return;
    }
    LFatal("This function can not be used for Array<realArray>.");
}

template <>
void OpenHurricane::cutProcessReceiv<OpenHurricane::realArray>::recvNonBlock(HurMPI::Request *request,
                                                                     const integer sizeBuf) {
    buf_.resize(sizeBuf, Zero);
    auto ierr = HurMPI::irecv(buf_.data(), buf_.size(), feature<real>::MPIType, zone_.sendProc(),
                              zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Receiv data error: %s", estr.c_str());
    }
}

template <>
void OpenHurricane::perProcessReceiv<OpenHurricane::realArray>::recvNonBlock(HurMPI::Request *request,
                                                                     const integer sizeBuf) {
    if (zone_.sendProc() == zone_.receivProc()) {
        return;
    }

    buf_.resize(sizeBuf, Zero);
    auto ierr = HurMPI::irecv(buf_.data(), buf_.size(), feature<real>::MPIType, zone_.sendProc(),
                              zone_.sendProc(), HurMPI::getComm(), request);

    if (ierr != MPI_SUCCESS) {
        std::string estr;
        HurMPI::errorString(ierr, estr);
        LFatal("Receiv data error: %s", estr.c_str());
    }
}

//==================================================================================
// For vector

template <> void OpenHurricane::perProcessReceiv<OpenHurricane::vector>::readBUf(Array<vector> &f) {
    if (buf_.size() != feature<vector>::nElements_ * zone_.des().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        for (int j = 0; j < feature<vector>::nElements_; j++) {
            f[zone_.des()[i]][j] = buf_[i * feature<vector>::nElements_ + j];
        }

        zone_.rotate(f[zone_.des()[i]]);
    }
}

//==================================================================================
// For tensor

template <> void OpenHurricane::perProcessReceiv<OpenHurricane::tensor>::readBUf(Array<tensor> &f) {
    if (buf_.size() != feature<tensor>::nElements_ * zone_.des().size()) {
        LFatal("The size of buf is less than the size of des.");
    }

    for (integer i = 0; i < zone_.des().size(); i++) {
        for (int j = 0; j < feature<tensor>::nElements_; j++) {
            f[zone_.des()[i]][j] = buf_[i * feature<tensor>::nElements_ + j];
        }

        zone_.rotate(f[zone_.des()[i]]);
    }
}
