/*!
 * \file thermo.cpp
 * \brief Main subroutines for thermo properties.
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
#include "thermo.hpp"

const OpenHurricane::real OpenHurricane::thermo::tol_ = 1.0e-4;
const int OpenHurricane::thermo::maxIter_ = 100;

namespace OpenHurricane {
    createClassNameStr(thermo, "thermo");
}

namespace OpenHurricane {
    createObjFty(thermo, controller);
}

OpenHurricane::thermo::thermo(const controller &cont, const equationOfState &EOS, const integer id)
    : eosPtr_(&EOS), phase_(phaseType::gas), speciesIndex_(id) {
    string pType;
    if (cont.found("phaseType", true)) {
        pType = cont.findWord("phaseType", true);
        stringToUpperCase(pType);
    }

    if (pType == "G") {
        phase_ = phaseType::gas;
    } else if (pType == "S") {
        phase_ = phaseType::solid;
    } else if (pType == "L") {
        phase_ = phaseType::liquid;
    }
}

OpenHurricane::thermo::thermo(const equationOfState &EOS, const integer id)
    : eosPtr_(&EOS), phase_(phaseType::gas), speciesIndex_(id) {}

OpenHurricane::thermo::thermo(const equationOfState &EOS, const integer id, const phaseType pt)
    : eosPtr_(&EOS), phase_(pt), speciesIndex_(id) {}

OpenHurricane::thermo::thermo(const equationOfState &EOS, const integer id, const std::string &pt)
    : eosPtr_(&EOS), phase_(phaseType::gas), speciesIndex_(id) {
    std::string pType = stringToUpperCase(pt);
    if (pType == "G") {
        phase_ = phaseType::gas;
    } else if (pType == "S") {
        phase_ = phaseType::solid;
    } else if (pType == "L") {
        phase_ = phaseType::liquid;
    }
}

OpenHurricane::thermo::thermo(const thermo &t)
    : eosPtr_(t.eosPtr_), phase_(t.phase_), speciesIndex_(t.speciesIndex_) {}


OpenHurricane::thermo &OpenHurricane::thermo::operator=(const thermo &t) {
    if (this != std::addressof(t)) {
        speciesIndex_ = t.speciesIndex_;
        phase_ = t.phase_;
    }
    return *this;
}

OpenHurricane::uniquePtr<OpenHurricane::thermo>
OpenHurricane::thermo::creator(const controller &cont, const equationOfState &EOS, const integer id) {
    string thermoType = cont.findWord("type");

    if (cont.found(thermoType)) {
        defineInObjCreator(thermo, thermoType, controller,
                           (cont.subController(thermoType), EOS, id));
    } else {
        defineInObjCreator(thermo, thermoType, controller, (cont, EOS, id));
    }
}

OpenHurricane::real OpenHurricane::thermo::T(real f, real p, real T0, integer &iFlag,
                                     real (thermo::*F)(const real, const real) const,
                                     real (thermo::*dFdT)(const real, const real) const,
                                     real (thermo::*limit)(const real, integer &) const) const {

    real Test = T0;
    real Tnew = T0;
    real Ttol = T0 * tol_;
    int iter = 0;

    do {
        Test = Tnew;
        Tnew = (this->*limit)(Test - ((this->*F)(p, Test) - f) / (this->*dFdT)(p, Test), iFlag);

        if (iter++ > maxIter_) {
#ifdef HUR_DEBUG
            PLWarning("Maximum number of iterations exceeded: %d", maxIter_);
#else
            if (report) {
                LWarning("Maximum number of iterations exceeded: %d", maxIter_);
            }
#endif // HUR_DEBUG

            break;
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}
