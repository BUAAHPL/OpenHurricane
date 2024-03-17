/*!
 * \file Pout.hpp
 * \brief Header of program runtime information
 *       The subroutines and functions are in the <i>Pout.cpp</i> file.
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

#include "HurMPIBase.hpp"
#include "fileName.hpp"
#include "real.hpp"
#include "string.hpp"
#include <iostream>
#include "HurFormat.hpp"

namespace OpenHurricane {

    /**
     * \brief The class of printing information to the screen.
     */
    class printToScreen {
    private:
        std::streamsize defaultCoutPrecision_;

    public:
        printToScreen() : defaultCoutPrecision_(std::cout.precision()) {}

        void greeting(const int proc = HurMPIBase::masterNo());

        void printHelp(const int proc = HurMPIBase::masterNo());

        template <typename T> inline printToScreen &operator<<(const T &);

        inline printToScreen &operator<<(std::ostream &(*op)(std::ostream &));

        inline printToScreen &operator<<(printToScreen &);

        printToScreen &setReal(const int _precision = feature<real>::precision);
        printToScreen &unsetReal();

        printToScreen &setScientic(const int _precision = feature<real>::precision);
        printToScreen &unsetScientic();

        int setf(int flag) const;
        void unsetf(int flag) const;

        std::streamsize width() const;
        std::streamsize width(std::streamsize w) const;

        printToScreen &fill(char c);
        printToScreen &fillDefault();

        template <typename... Args> inline void out(const char* fmt, Args... args) {
            if (HurMPIBase::master()) {
                std::cout << hurFormat(fmt, args...);
            }
        }

        template <typename... Args> inline void operator()(const char *fmt, Args... args) {
            out(fmt, args...);
        }
    };

    /**
     * \brief Print to screen on master process.
     */
    static printToScreen Pout;

} // namespace OpenHurricane

#include "Pout.inl"