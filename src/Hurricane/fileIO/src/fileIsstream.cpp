/*!
 * \file fileIsstream.cpp
 * \brief Main subroutines of file input string stream.
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

#include "fileIsstream.hpp"
#include "dataStructure.hpp"

OpenHurricane::fileIsstream::fileIsstream() : IOsstream(), bufferPtr_(nullptr), isBuffered_(false) {}

OpenHurricane::fileIsstream::fileIsstream(const streamFormat format, std::ios_base::openmode mode)
    : IOsstream(format, mode), bufferPtr_(nullptr), isBuffered_(false) {}

OpenHurricane::fileIsstream::fileIsstream(const fileName &name, const streamFormat format,
                                      std::ios_base::openmode mode)
    : IOsstream(name, format, mode), bufferPtr_(nullptr), isBuffered_(false) {
    open();
}

void OpenHurricane::fileIsstream::open() {
    if (opened()) {
        if (isBuffered_) {
            LFatal("Attempt to open an opened file in this stream. File name: %s", name().c_str());
        }
    }

    if (!name().empty()) {
        if (!existed()) {
            LFatal("Attempt to open an inexistent file: %s", name().c_str());
        }

        if_.open(name().c_str(), mode());
        setOpened();
    } else {
        LFatal("Attempt to open a null name file.");
    }
}

void OpenHurricane::fileIsstream::close() {
    if (opened()) {
        if_.close();
        setClosed();
    }
}

void OpenHurricane::fileIsstream::readBuffer() {
    if (isBuffered_) {
        LFatal("Attempt to buffer the file: %s again.", name().c_str());
    }

    if (closed()) {
        open();
    }

    clear();

    char *buffer;
    std::streampos lenth = if_.seekg(0, std::ios::end).tellg();
    buffer = new char[static_cast<std::streamsize>(lenth) + 1];

    if_.seekg(0, std::ios::beg).read(buffer, static_cast<std::streamsize>(lenth));

    buffer[static_cast<std::streamsize>(lenth)] = '\0';

    bufferPtr_.reset(new std::string(buffer));

    delete[] buffer;

    isBuffered_ = true;
}