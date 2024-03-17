/*!
 * \file compactSupportFunction.cpp
 * \brief The subroutines and functions of compact support function for RBF.
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

#include "compactSupportFunction.hpp"

namespace OpenHurricane {
    const std::string compactSupportFunction::className_ = "compactSupportFunction";
    createObjFty(compactSupportFunction, controller);

    compactSupportFunction::compactSupportFunction(const controller &cont)
        : R_(cont.findOrDefault<real>("R", 1.0)) {}

    uniquePtr<compactSupportFunction> compactSupportFunction::creator(const controller &cont) {
        string compactSupportFunctionType = cont.findWord(compactSupportFunction::className_);
        Pout << "    Setting compact support function: " << compactSupportFunctionType << std::endl;
        defineInObjCreator(compactSupportFunction, compactSupportFunctionType, controller, (cont));
    }
} // namespace OpenHurricane

namespace OpenHurricane {
    const std::string CPC0::className_ = "CPC0";
    registerObjFty(compactSupportFunction, CPC0, controller);

    CPC0::CPC0(const controller &cont) : compactSupportFunction(cont) {}
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::CPC0::f(const real psi) const {
    if (psi >= 1.0) {
        return 0.0;
    } else {
        return sqr(real(1) - psi);
    }
}

namespace OpenHurricane {
    const std::string CPC2::className_ = "CPC2";
    registerObjFty(compactSupportFunction, CPC2, controller);

    CPC2::CPC2(const controller &cont) : compactSupportFunction(cont) {}
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::CPC2::f(const real psi) const {
    if (psi >= 1.0) {
        return 0.0;
    } else {
        return pow4(real(1.0) - psi) * (4.0 * psi + 1.0);
    }
}

namespace OpenHurricane {
    const std::string CPC4::className_ = "CPC4";
    registerObjFty(compactSupportFunction, CPC4, controller);

    CPC4::CPC4(const controller &cont) : compactSupportFunction(cont) {}
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::CPC4::f(const real psi) const {
    if (psi >= 1.0) {
        return 0.0;
    } else {
        return pow6(real(1) - psi) * (35.0 / 3.0 * sqr(psi) + 6.0 * psi + 1.0);
    }
}

namespace OpenHurricane {
    const std::string CPC6::className_ = "CPC6";
    registerObjFty(compactSupportFunction, CPC6, controller);

    CPC6::CPC6(const controller &cont) : compactSupportFunction(cont) {}
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::CPC6::f(const real psi) const {
    if (psi >= 1.0) {
        return 0.0;
    } else {
        return sqr(pow4(real(1) - psi)) * (32.0 * pow3(psi) + 25.0 * sqr(psi) + 8.0 * psi + 1.0);
    }
}

namespace OpenHurricane {
    const std::string CTPSC0::className_ = "CTPSC0";
    registerObjFty(compactSupportFunction, CTPSC0, controller);

    CTPSC0::CTPSC0(const controller &cont) : compactSupportFunction(cont) {}
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::CTPSC0::f(const real psi) const {
    if (psi >= 1.0) {
        return 0.0;
    } else {
        return pow5(real(1) - psi);
    }
}

namespace OpenHurricane {
    const std::string CTPSC1::className_ = "CTPSC1";
    registerObjFty(compactSupportFunction, CTPSC1, controller);

    CTPSC1::CTPSC1(const controller &cont) : compactSupportFunction(cont) {}
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::CTPSC1::f(const real psi) const {
    if (psi >= 1.0) {
        return 0.0;
    } else {
        return 1.0 + 80.0 / 3.0 * sqr(psi) - 40.0 * pow3(psi) + 15.0 * pow4(psi) -
               8.0 / 3.0 * pow5(psi) + 20.0 * sqr(psi) * log(psi);
    }
}

namespace OpenHurricane {
    const std::string CTPSC2a::className_ = "CTPSC2a";
    registerObjFty(compactSupportFunction, CTPSC2a, controller);

    CTPSC2a::CTPSC2a(const controller &cont) : compactSupportFunction(cont) {}
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::CTPSC2a::f(const real psi) const {
    if (psi >= 1.0) {
        return 0.0;
    } else {
        return 1.0 - 30.0 * sqr(psi) - 10.0 * pow3(psi) + 45.0 * pow4(psi) - 6.0 * pow5(psi) -
               60.0 * pow3(psi) * log(psi);
    }
}

namespace OpenHurricane {
    const std::string CTPSC2b::className_ = "CTPSC2b";
    registerObjFty(compactSupportFunction, CTPSC2b, controller);

    CTPSC2b::CTPSC2b(const controller &cont) : compactSupportFunction(cont) {}
} // namespace OpenHurricane

hur_nodiscard OpenHurricane::real OpenHurricane::CTPSC2b::f(const real psi) const {
    if (psi >= 1.0) {
        return 0.0;
    } else {
        return 1.0 - 20.0 * sqr(psi) + 80.0 * pow3(psi) - 45.0 * pow4(psi) - 16.0 * pow5(psi) +
               60.0 * pow4(psi) * log(psi);
    }
}