/*!
 * \file fileOsstream.cpp
 * \brief Main subroutines of file output string stream.
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

#include "fileOsstream.hpp"
#include "dataStructure.hpp"
#include "logFile.hpp"

OpenHurricane::fileOsstream::fileOsstream()
    : IOsstream(), defaultCoutPrecision_(of_.precision()), openOption_(CREATE_MASTER) {}

OpenHurricane::fileOsstream::fileOsstream(const streamFormat format, std::ios_base::openmode mode)
    : IOsstream(format, mode), defaultCoutPrecision_(of_.precision()), openOption_(CREATE_MASTER) {}

OpenHurricane::fileOsstream::fileOsstream(const streamFormat format, const openOptions _op,
                                      std::ios_base::openmode mode)
    : IOsstream(format, mode), defaultCoutPrecision_(of_.precision()), openOption_(_op) {}

OpenHurricane::fileOsstream::fileOsstream(const fileName &name, const streamFormat format,
                                      std::ios_base::openmode mode)
    : IOsstream(name, format, mode), defaultCoutPrecision_(of_.precision()),
      openOption_(CREATE_MASTER) {
    open();
}

OpenHurricane::fileOsstream::fileOsstream(const fileName &name, const openOptions _op,
                                      const streamFormat format, std::ios_base::openmode mode)
    : IOsstream(name, format, mode), defaultCoutPrecision_(of_.precision()), openOption_(_op) {
    open();
}

void OpenHurricane::fileOsstream::open() {
    if (opened()) {
        LFatal("Attempt to open an opened file in this stream. File name: %s", name().c_str());
    }

    if (!name().empty()) {
        if (openOption_ == CREATE_ALL) {
            of_.open(name().c_str(), mode());
            if (!of_.fail()) {
                setOpened();
                Pout("    Info: Creatting file: \"%s\" in all processors.\n", name().c_str());
                if (report) {
                    LInfo("Creating file: \"%s\" in all processors", name().c_str());
                }
            }
        } else if (openOption_ == CREATE_MASTER) {
            if (HurMPIBase::master()) {
                of_.open(name().c_str(), mode());
                if (!of_.fail()) {
                    setOpened();
                    if (report) {
                        LInfo("Creating file: \"%s\" in all master", name().c_str());
                    }
                }
            }
            HurMPIBase::barrier(HurMPIBase::getComm());
            if (!HurMPIBase::master()) {
                std::fstream _file;
                //_file.open(name().c_str(), std::ios::in | std::ios::_Nocreate);
                _file.open(name().c_str(), std::ios::in);
                if (!_file.fail()) {
                    setOpened();
                }
                _file.close();

                if (opened()) {
                    of_.open(name().c_str(), mode() | std::ios::app);
                    if (report) {
                        LInfo("Openning file: \"%s\" in slavers", name().c_str());
                    }
                } else {
                    setClosed();
                }
            }
        } else {
            if (HurMPIBase::master()) {
                of_.open(name().c_str(), mode());
                if (!of_.fail()) {
                    setOpened();
                    if (report) {
                        LInfo("Creating file: \"%s\" only in slavers", name().c_str());
                    }
                }
            }
        }
    } else {
        LFatal("Attempt to open a null name file.");
    }
}

void OpenHurricane::fileOsstream::close() {
    if (opened()) {
        of_.close();
        setClosed();
    }
}

std::ofstream &OpenHurricane::fileOsstream::flush() {
    if (opened()) {
        of_.flush();
    }
    return of_;
}

hur_nodiscard std::ofstream &OpenHurricane::fileOsstream::os() noexcept {
    return of_;
}

std::ofstream &OpenHurricane::fileOsstream::write(const char *_str, std::streamsize _count) {
    of_.write(_str, _count);
    return of_;
}

std::ofstream &OpenHurricane::fileOsstream::write(const std::string &_str, std::streamsize _count) {
    of_.write(_str.c_str(), _count);
    return of_;
}

std::ofstream &OpenHurricane::fileOsstream::setRealPrecision(const int _precision) {
    of_.precision(_precision);
    return of_;
}

std::ofstream &OpenHurricane::fileOsstream::unsetRealPrecision() {
    of_.precision(defaultCoutPrecision_);
    return of_;
}

std::ofstream &OpenHurricane::fileOsstream::setScientic(const int _precision) {
    of_.setf(std::ios::scientific);
    of_.precision(_precision);
    return of_;
}

std::ofstream &OpenHurricane::fileOsstream::unsetScientic() {
    of_.unsetf(std::ios::scientific);
    of_.precision(defaultCoutPrecision_);
    return of_;
}
