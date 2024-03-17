#include "AUSMPlusUP.hpp"
/*!
 * \file AUSMPlusUP.inl
 * \brief The In-Line functions of the <i>AUSMPlusUP.hpp</i> file.
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
#include <cmath>

hur_nodiscard inline OpenHurricane::real OpenHurricane::AUSMPlusUP::alpha(const real fa) const {
    return 3.0 / 16.0 * (-4.0 + 5.0 * fa * fa);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::AUSMPlusUP::fa(const real Mo) const {
    return Mo * (2.0 - Mo);
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::AUSMPlusUP::M1P(const real Ma) const {
    return 0.5 * (Ma + mag(Ma));
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::AUSMPlusUP::M1M(const real Ma) const {
    return 0.5 * (Ma - mag(Ma));
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::AUSMPlusUP::M2P(const real Ma) const {
    return 0.25 * sqr(Ma + real(1.0));
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::AUSMPlusUP::M2M(const real Ma) const {
    return -0.25 * sqr(Ma - real(1.0));
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::AUSMPlusUP::M4P(const real Ma) const {
    real m;
    if (mag(Ma) >= 1.0) {
        m = M1P(Ma);
    } else {
        m = M2P(Ma) * (1.0 - 16.0 * beta_ * M2M(Ma));
    }
    return m;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::AUSMPlusUP::M4M(const real Ma) const {
    real m;
    if (mag(Ma) >= 1.0) {
        m = M1M(Ma);
    } else {
        m = M2M(Ma) * (1.0 + 16.0 * beta_ * M2P(Ma));
    }
    return m;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::AUSMPlusUP::P5P(const real Ma,
                                                                        const real _alpha) const {
    real p;
    if (mag(Ma) >= 1.0) {
        p = M1P(Ma) / Ma;
    } else {
        p = M2P(Ma) * ((2.0 - Ma) - 16.0 * _alpha * Ma * M2M(Ma));
    }
    return p;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::AUSMPlusUP::P5M(const real Ma,
                                                                        const real _alpha) const {
    real p;
    if (mag(Ma) >= 1.0) {
        p = M1M(Ma) / Ma;
    } else {
        p = M2M(Ma) * ((-2.0 - Ma) + 16.0 * _alpha * Ma * M2P(Ma));
    }
    return p;
}
