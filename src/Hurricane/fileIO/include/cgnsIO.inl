#include "cgnsIO.hpp"
/*!
 * \file cgnsIO.inl
 * \brief In-Line subroutines of the <i>cgnsIO.hpp</i> file.
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
#ifdef USES_CGNS

inline void OpenHurricane::cgnsIO::setOpened() noexcept {
    openClosed_ = OPENED;
}

inline void OpenHurricane::cgnsIO::setClosed() noexcept {
    openClosed_ = CLOSED;
}

inline OpenHurricane::cgnsIO::cgnsIO() : filename_(), cgnsFN_(-1), mode_(MODIFY), openClosed_(CLOSED) {}

inline OpenHurricane::cgnsIO::cgnsIO(const fileName &fN)
    : filename_(fN), cgnsFN_(-1), mode_(MODIFY), openClosed_(CLOSED) {}

inline OpenHurricane::cgnsIO::~cgnsIO() noexcept {
    if (opened()) {
        close();
    }
}

inline void OpenHurricane::cgnsIO::open(const fileName &fN) {
    filename_ = fN;
    open();
}

inline void OpenHurricane::cgnsIO::open(const unsigned int flg) {
    setMode(flg);
    open();
}

hur_nodiscard inline float OpenHurricane::cgnsIO::cgVersion() {
    float ver = 0.0f;
    const auto result = cg_version(cgnsFN_, &ver);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return ver;
}

hur_nodiscard inline int OpenHurricane::cgnsIO::cgPrecision() {
    int ip = 0;
    const auto result = cg_precision(cgnsFN_, &ip);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return ip;
}

hur_nodiscard inline int OpenHurricane::cgnsIO::isCGNS() {
    int ip = 0;
    const auto result = cg_is_cgns(filename_.c_str(), &ip);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return ip;
}

inline void OpenHurricane::cgnsIO::setCGNSFileType(const int ft) {
    const auto result = cg_set_file_type(ft);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
}

hur_nodiscard inline int OpenHurricane::cgnsIO::getCGNSFileType() {
    int ip = 0;
    const auto result = cg_get_file_type(cgnsFN_, &ip);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return ip;
}

inline void OpenHurricane::cgnsIO::close() {
    if (opened()) {
        cg_close(cgnsFN_);
        setClosed();
    }
}

hur_nodiscard inline bool OpenHurricane::cgnsIO::opened() const noexcept {
    return openClosed_ == OPENED;
}

hur_nodiscard inline bool OpenHurricane::cgnsIO::closed() const noexcept {
    return openClosed_ == CLOSED;
}

hur_nodiscard inline bool OpenHurricane::cgnsIO::isRead() const noexcept {
    return mode_ == READ;
}

hur_nodiscard inline bool OpenHurricane::cgnsIO::isWrite() const noexcept {
    return mode_ == WRITE;
}

hur_nodiscard inline bool OpenHurricane::cgnsIO::isModify() const noexcept {
    return mode_ == MODIFY;
}

hur_nodiscard inline DataType_t OpenHurricane::cgnsIO::getDataTypeFromTheProgram() const {
    if (feature<real>::dataFormat == 1) {
        return RealSingle;
    } else if (feature<real>::dataFormat == 2) {
        return RealDouble;
    } else {
        LFatal("Unsupported data format");
    }
    return RealDouble;
}

hur_nodiscard inline int OpenHurricane::cgnsIO::readCellDim(const int B) const {
    int cellDim;
    readCellDim(B, cellDim);
    return cellDim;
}

inline int OpenHurricane::cgnsIO::writeBase(const string &baseName) {
    return writeBase(baseName, 3, 3);
}

#endif // USES_CGNS