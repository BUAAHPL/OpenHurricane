#include "IOsstream.hpp"
/*!
 * \file IOsstream.inl
 * \brief In-Line subroutines of the <i>IOsstream.hpp</i> file.
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

inline void OpenHurricane::IOsstream::setOpened() { openClosed_ = OPENED; }

inline void OpenHurricane::IOsstream::setClosed() { openClosed_ = CLOSED; }

inline OpenHurricane::IOsstream::IOsstream()
    : name_(), format_(ASCII_FORMAT), mode_(std::ios_base::in),
      openClosed_(CLOSED) {}

inline OpenHurricane::IOsstream::IOsstream(const streamFormat format,
                                           std::ios_base::openmode mode)
    : name_(), format_(format), mode_(mode), openClosed_(CLOSED){}

inline OpenHurricane::IOsstream::IOsstream(const fileName &name,
                                           const streamFormat format,
                                           std::ios_base::openmode mode)
    : name_(name), format_(format), mode_(mode), openClosed_(CLOSED) {}

hur_nodiscard inline const OpenHurricane::fileName &OpenHurricane::IOsstream::name() const noexcept {
    return name_;
}

hur_nodiscard inline OpenHurricane::fileName &OpenHurricane::IOsstream::name() noexcept {
    return name_;
}

hur_nodiscard inline OpenHurricane::IOsstream::streamFormat
OpenHurricane::IOsstream::format() const noexcept {
    return format_;
}

hur_nodiscard inline std::ios_base::openmode OpenHurricane::IOsstream::mode() const noexcept {
    if (format_ == ASCII_FORMAT) {
        return mode_;
    } else {
        return mode_ | std::ios_base::binary;
    }
}

hur_nodiscard inline bool OpenHurricane::IOsstream::opened() const noexcept {
    return openClosed_ == OPENED;
}

hur_nodiscard inline bool OpenHurricane::IOsstream::closed() const noexcept {
    return openClosed_ == CLOSED;
}

inline void OpenHurricane::IOsstream::changeMode(const std::ios_base::openmode _mode) noexcept {
    mode_ = _mode;
}

inline void OpenHurricane::IOsstream::changeFileName(const fileName &_name) noexcept {
    name_ = _name;
}

inline void OpenHurricane::IOsstream::changeFormat(const streamFormat _format) noexcept {
    format_ = _format;
}

hur_nodiscard inline bool OpenHurricane::IOsstream::existed() const {
    std::fstream _file;
    //_file.open(name_.c_str(), std::ios::in | std::ios::_Nocreate);
    _file.open(name_.c_str(), std::ios::in);
    if (_file.fail()) {
        _file.close();
        return false;
    }
    _file.close();
    return true;
}
