/*!
 * \file Pout.cpp
 * \brief The subroutines and functions of program information
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

#include "Pout.hpp"
#include "GPUInfo.hpp"
#include "HurMPIBase.hpp"
#include "clock.hpp"
#include "version.hpp"
#include <iomanip>

void OpenHurricane::printToScreen::greeting(const int proc) {
    if (HurMPIBase::isThisProc(proc)) {
        std::string timeName = "     Starting time: ";
        timeName += clock::dateTime();
        std::cout << "    =" << std::setw(65) << std::setfill('=') << "="
                  << std::endl;
        std::cout
            << std::left << "    ||" << std::setfill(' ') << std::setw(62)
            << "  OpenHurricane: Unstructured Computational Fluid Dynamics Program"
            << "||" << std::endl;
        std::cout << "    ||" << std::setw(62)
                  << "     Copyright(c) 2019-2022 BUAA XX Research Group"
                  << "||" << std::endl;
        std::cout << "    ||" << std::setw(62) << timeName.c_str() << "||"
                  << std::endl;
        std::cout << "    =" << std::setw(65) << std::setfill('=') << "="
                  << std::endl;
        std::cout << std::right << std::setfill(' ');
    }
    HurMPIBase::barrier(HurMPIBase::getComm());
    HurMPIBase::getPcInfo();
    GPUInfo::getGPUInfo();
}

void OpenHurricane::printToScreen::printHelp(const int proc) {}

OpenHurricane::printToScreen &
OpenHurricane::printToScreen::setReal(const int _precision) {
    std::cout.precision(_precision);
    return *this;
}

OpenHurricane::printToScreen &OpenHurricane::printToScreen::unsetReal() {
    std::cout.precision(defaultCoutPrecision_);
    return *this;
}

OpenHurricane::printToScreen &
OpenHurricane::printToScreen::setScientic(const int _precision) {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(_precision);
    return *this;
}

OpenHurricane::printToScreen &OpenHurricane::printToScreen::unsetScientic() {
    std::cout.unsetf(std::ios::scientific);
    std::cout.precision(defaultCoutPrecision_);
    return *this;
}

int OpenHurricane::printToScreen::setf(int flag) const {
    return int(std::cout.setf(std::ios_base::fmtflags(flag)));
}

void OpenHurricane::printToScreen::unsetf(int flag) const {
    std::cout.unsetf(std::ios_base::fmtflags(flag));
}

OpenHurricane::printToScreen &OpenHurricane::printToScreen::fill(char c) {
    std::cout << std::setfill(c);
    return *this;
}

OpenHurricane::printToScreen &OpenHurricane::printToScreen::fillDefault() {
    std::cout << std::setfill(' ');
    return *this;
}
