/*!
 * \file testProcessTime.cpp
 * \brief Main subroutines for testing the process time of OpenHurricane.
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
#include "testProcessTime.hpp"
#include "iteration.hpp"

OpenHurricane::testProcessTime::testProcessTime(const iteration &iter, const controller &cont)
    : fosTestPT_(IOsstream::streamFormat::ASCII_FORMAT, IOsstream::openOptions::ONLY_MASTER,
                 std::ios_base::out),
      hClockPtr_(nullptr), titleBuffStr_(), firstLineBuffStr_(), count_(0) {
    if (cont.found("fileName")) {
        string pw = cont.findWord("fileName");
        trim(pw);
        fileName outFile = pw;
        if (!outFile.isAbsolute()) {
            outFile = iter.outputPath() / outFile;
        }

        fileName resOut = outFile;
        fosTestPT_.changeFileName(resOut);

        fosTestPT_.open();
    } else {
        LFatal("Cannot fing fileName in %s", cont.name().c_str());
    }
}

OpenHurricane::testProcessTime::testProcessTime(const iteration &iter, const fileName &fn)
    : fosTestPT_(IOsstream::streamFormat::ASCII_FORMAT, IOsstream::openOptions::ONLY_MASTER,
                 std::ios_base::out),
      hClockPtr_(nullptr), titleBuffStr_(), firstLineBuffStr_(), count_(0) {
    fileName outFile = fn;
    if (!outFile.isAbsolute()) {
        outFile = iter.outputPath() / outFile;
    }

    fileName resOut = outFile;
    fosTestPT_.changeFileName(resOut);

    fosTestPT_.open();
}

void OpenHurricane::testProcessTime::start(const integer istep) const {
    if (count_ > 0) {
        fosTestPT_.os() << istep;
    } else {
        titleBuffStr_.clear();
        titleBuffStr_ = "\"step\",";
        firstLineBuffStr_.clear();
        firstLineBuffStr_ = toString(istep);
        firstLineBuffStr_ += "\t";
    }
    HurMPIBase::barrier();
    start();
}

hur_nodiscard OpenHurricane::real OpenHurricane::testProcessTime::clockTimeIncrement() const {
    HurMPIBase::barrier();
    if (hClockPtr_) {
        return hClockPtr_->clockTimeIncrement();
    }
    return static_cast<real>(0);
}

void OpenHurricane::testProcessTime::clockTimeIncrement(const char *name) const {
    HurMPIBase::barrier();
    if (hClockPtr_) {
        if (count_ > 0) {
            fosTestPT_.os() << '\t' << std::setprecision(7) << hClockPtr_->clockTimeIncrement();
        } else {
            titleBuffStr_ += "\"";
            titleBuffStr_ += name;
            titleBuffStr_ += "\",";
            firstLineBuffStr_ += toString(hClockPtr_->clockTimeIncrement());
            firstLineBuffStr_ += "\t";
        }
    }
}

void OpenHurricane::testProcessTime::clockTimeIncrement(const char *name, const real time) const {
    HurMPIBase::barrier();
    if (hClockPtr_) {
        if (count_ > 0) {
            fosTestPT_.os() << '\t' << std::setprecision(7) << time;
        } else {
            titleBuffStr_ += "\"";
            titleBuffStr_ += name;
            titleBuffStr_ += "\",";
            firstLineBuffStr_ += toString(time);
            firstLineBuffStr_ += "\t";
        }
    }
}

void OpenHurricane::testProcessTime::stop() const {
    HurMPIBase::barrier();
    if (hClockPtr_) {
        if (count_ > 0) {
            fosTestPT_.os() << '\t' << std::setprecision(7) << hClockPtr_->elapsedClockTime()
                            << std::endl;
        } else {
            fosTestPT_.os() << "TITLE=\"OpenHurricane testing process time\"" << std::endl;
            fosTestPT_.os() << "variables = ";
            fosTestPT_.os() << titleBuffStr_ << "\"elapsed\"" << std::endl;
            fosTestPT_.os() << firstLineBuffStr_ << hClockPtr_->elapsedClockTime() << std::endl;
            titleBuffStr_.clear();
            firstLineBuffStr_.clear();
        }
    }
    resetHZClock();
    count_++;
}