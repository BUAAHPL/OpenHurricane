/*!
 * \file makeDualTimeSteppingMethod.cpp
 * \brief Main subroutines for making dual time-stepping method.
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
#include "makeDualTimeSteppingMethod.hpp"
#include "solver.hpp"

// Dual time-stepping with LU-SGS

namespace OpenHurricane {
    createClassNameStrTmpl(BDF123LUSGS, "BDF123LUSGS");
    createClassNameStrTmpl(ESDIRKLUSGS, "ESDIRKLUSGS");
    registerObjFty(timeMarching, BDF123LUSGS, controller);
    registerObjFty(timeMarching, ESDIRKLUSGS, controller);
} // namespace OpenHurricane

template <> void OpenHurricane::BDF123<OpenHurricane::LUSGS>::setSolverWrite() {
    timeMarching::solver_.setBDFUnsteadySolver();
    const_cast<iteration &>(timeMarching::iter()).setWriteLastToRelay(true);
    const_cast<iteration &>(timeMarching::iter()).setReadLastFromRelay(true);
}

template <>
void OpenHurricane::BDF123<OpenHurricane::LUSGS>::setBDFType(const controller &timeMethodCont) {
    if (timeMethodCont.subController("BDF123LUSGS").found("BDF123")) {
        string t = timeMethodCont.subController("BDF123LUSGS").findWord("BDF123");

        if (t == "BDF1") {
            type_ = BDF1;
        } else if (t == "BDF2") {
            type_ = BDF2;
        } else if (t == "BDF3") {
            type_ = BDF3;
        } else {
            errorAbortStr(("Unknown BDF type: " + t));
        }
    }
}

template <> void OpenHurricane::ESDIRK<OpenHurricane::LUSGS>::setSolverWrite() {
    const_cast<iteration &>(timeMarching::iter()).setWriteLastToRelay(true);
    const_cast<iteration &>(timeMarching::iter()).setReadLastFromRelay(true);
}

template <>
void OpenHurricane::ESDIRK<OpenHurricane::LUSGS>::setESDIRKType(const controller &timeMethodCont) {
    if (timeMethodCont.subController("ESDIRKLUSGS").found("ESDIRK")) {
        string t = timeMethodCont.subController("ESDIRKLUSGS").findWord("ESDIRK");

        // Four stages third order ESDIRK
        if (t == "ESDIRK3") {
            type_ = ESDIRK3;
            stage_ = 4; //stage = 4
        }
        // Six stages fourth order ESDIRK
        else if (t == "ESDIRK4") {
            stage_ = 6; //stage = 6
            type_ = ESDIRK4;
        } else {
            errorAbortStr(("Unknown ESDIRK type: " + t));
        }
    }
}
