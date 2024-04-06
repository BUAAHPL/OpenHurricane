/*!
 * \file tensors.cpp
 * \brief Main subroutines for tensors.
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
#include "tensors.hpp"
#include "basicFunctions.hpp"
#include "mathConstants.hpp"

using namespace OpenHurricane::constant;

template <>
const std::streamsize OpenHurricane::tensor::vsType::precision(feature<real>::precision);

template <> const OpenHurricane::tensor OpenHurricane::tensor::I(1, 0, 0, 0, 1, 0, 0, 0, 1);

template <>
const std::streamsize OpenHurricane::integerTensor::vsType::precision(feature<integer>::precision);

#ifdef HURRICANE_DP
#include "floatReal.hpp"
template <>
const std::streamsize OpenHurricane::floatTensor::vsType::precision(feature<floatReal>::precision);
template <>
const OpenHurricane::floatTensor OpenHurricane::floatTensor::I(1, 0, 0, 0, 1, 0, 0, 0, 1);
#endif //HURRICANE_DP

template <>
const std::streamsize OpenHurricane::symmTensor::vsType::precision(feature<real>::precision);
template <> const OpenHurricane::symmTensor OpenHurricane::symmTensor::I(1, 0, 0, 1, 0, 1);

template <>
const std::streamsize
    OpenHurricane::integerSymmTensor::vsType::precision(feature<integer>::precision);
template <>
const OpenHurricane::integerSymmTensor OpenHurricane::integerSymmTensor::I(1, 0, 0, 1, 0, 1);

template <>
const std::streamsize OpenHurricane::diagTensor::vsType::precision(feature<real>::precision);
template <> const OpenHurricane::diagTensor OpenHurricane::diagTensor::I(1, 1, 1);

template <>
const std::streamsize
    OpenHurricane::integerDiagTensor::vsType::precision(feature<integer>::precision);
template <> const OpenHurricane::integerDiagTensor OpenHurricane::integerDiagTensor::I(1, 1, 1);

template <>
const std::streamsize OpenHurricane::sphericalTensor::precision(feature<real>::precision);

template <> const OpenHurricane::sphericalTensor OpenHurricane::sphericalTensor::I(real(1));

template <>
const OpenHurricane::sphericalTensor OpenHurricane::sphericalTensor::oneThirdI(1.0 / 3.0);

template <>
const OpenHurricane::sphericalTensor OpenHurricane::sphericalTensor::twoThirdsI(2.0 / 3.0);

template <>
const std::streamsize OpenHurricane::integerSphericalTensor::precision(feature<integer>::precision);
template <> const OpenHurricane::integerSphericalTensor OpenHurricane::integerSphericalTensor::I(1);

template <> std::string OpenHurricane::toString(const tensor &st) {
    std::stringstream sxx;
    std::stringstream sxy;
    std::stringstream sxz;
    std::stringstream syx;
    std::stringstream syy;
    std::stringstream syz;
    std::stringstream szx;
    std::stringstream szy;
    std::stringstream szz;
    sxx << std::setprecision(feature<real>::precision) << st.xx();
    sxy << std::setprecision(feature<real>::precision) << st.xy();
    sxz << std::setprecision(feature<real>::precision) << st.xz();
    syx << std::setprecision(feature<real>::precision) << st.yx();
    syy << std::setprecision(feature<real>::precision) << st.yy();
    syz << std::setprecision(feature<real>::precision) << st.yz();
    szx << std::setprecision(feature<real>::precision) << st.zx();
    szy << std::setprecision(feature<real>::precision) << st.zy();
    szz << std::setprecision(feature<real>::precision) << st.zz();

    std::string s;
    s.append("\n");
    s += sxx.str() + ", " + sxy.str() + ", " + sxz.str();
    s.append("\n");
    s += syx.str() + ", " + syy.str() + ", " + syz.str();
    s.append("\n");
    s += szx.str() + ", " + szy.str() + ", " + szz.str();
    return s;
}

template <> hur_nodiscard std::string OpenHurricane::toString(const integerTensor &st) {
    std::string s;
    s.append("\n");
    s += std::to_string(st.xx()) + ", " + std::to_string(st.xy()) + ", " + std::to_string(st.xz());
    s.append("\n");
    s += std::to_string(st.yx()) + ", " + std::to_string(st.yy()) + ", " + std::to_string(st.yz());
    s.append("\n");
    s += std::to_string(st.zx()) + ", " + std::to_string(st.zy()) + ", " + std::to_string(st.zz());
    return s;
}

#ifdef HURRICANE_DP
template <> hur_nodiscard std::string OpenHurricane::toString(const floatTensor &st) {
    std::string s;
    s.append("\n");
    s += std::to_string(st.xx()) + ", " + std::to_string(st.xy()) + ", " + std::to_string(st.xz());
    s.append("\n");
    s += std::to_string(st.yx()) + ", " + std::to_string(st.yy()) + ", " + std::to_string(st.yz());
    s.append("\n");
    s += std::to_string(st.zx()) + ", " + std::to_string(st.zy()) + ", " + std::to_string(st.zz());
    return s;
}
#endif //HURRICANE_DP

template <> hur_nodiscard std::string OpenHurricane::toString(const symmTensor &st) {
    std::string s;
    s.append("\n");
    s += std::to_string(st.xx()) + ", " + std::to_string(st.xy()) + ", " + std::to_string(st.xz());
    s.append("\n");
    s += std::to_string(st.xy()) + ", " + std::to_string(st.yy()) + ", " + std::to_string(st.yz());
    s.append("\n");
    s += std::to_string(st.xz()) + ", " + std::to_string(st.yz()) + ", " + std::to_string(st.zz());
    return s;
}

template <> hur_nodiscard std::string OpenHurricane::toString(const integerSymmTensor &st) {
    std::string s;
    s.append("\n");
    s += std::to_string(st.xx()) + ", " + std::to_string(st.xy()) + ", " + std::to_string(st.xz());
    s.append("\n");
    s += std::to_string(st.xy()) + ", " + std::to_string(st.yy()) + ", " + std::to_string(st.yz());
    s.append("\n");
    s += std::to_string(st.xz()) + ", " + std::to_string(st.yz()) + ", " + std::to_string(st.zz());
    return s;
}

template <> hur_nodiscard std::string OpenHurricane::toString(const diagTensor &st) {
    std::string s;
    s.append("\n");
    s += std::to_string(st.xx()) + ", 0.0, 0.0";
    s.append("\n");
    s += "0.0, " + std::to_string(st.yy()) + ", 0.0";
    s.append("\n");
    s += "0.0, 0.0, " + std::to_string(st.zz());
    return s;
}

template <> hur_nodiscard std::string OpenHurricane::toString(const integerDiagTensor &st) {
    std::string s;
    s.append("\n");
    s += std::to_string(st.xx()) + ", 0.0, 0.0";
    s.append("\n");
    s += "0.0, " + std::to_string(st.yy()) + ", 0.0";
    s.append("\n");
    s += "0.0, 0.0, " + std::to_string(st.zz());
    return s;
}

template <> hur_nodiscard std::string OpenHurricane::toString(const sphericalTensor &st) {
    std::string s;
    s.append("\n");
    s += std::to_string(st.ii()) + ", 0.0, " + "0.0";
    s.append("\n");
    s += "0.0, " + std::to_string(st.ii()) + ", 0.0";
    s.append("\n");
    s += "0.0,  0.0, " + std::to_string(st.ii());
    return s;
}

template <> hur_nodiscard std::string OpenHurricane::toString(const integerSphericalTensor &st) {
    std::string s;
    s.append("\n");
    s += std::to_string(st.ii()) + ", 0, " + "0";
    s.append("\n");
    s += "0, " + std::to_string(st.ii()) + ", 0";
    s.append("\n");
    s += "0,  0, " + std::to_string(st.ii());
    return s;
}

hur_nodiscard OpenHurricane::real OpenHurricane::deltaij(const integer i, const integer j) {
#ifdef HUR_DEBUG
    if (i < 0 || i > 2 || j < 0 || j > 2) {
        errorAbortStr(
            ("The index (" + toString(i) + ", " + toString(j) + ") must be within (0~2, 0~2)"));
    }
#endif // HUR_DEBUG

    if (i == j) {
        return 1;
    }

    return 0;
}

hur_nodiscard OpenHurricane::real OpenHurricane::epsilonijk(const integer i, const integer j,
                                                            const integer k) {
#ifdef HUR_DEBUG
    if (i < 0 || i > 2 || j < 0 || j > 2 || k < 0 || k > 2) {
        errorAbortStr(("The index (" + toString(i) + ", " + toString(j) + ", " + toString(k) +
                       ") must be within (0~2, 0~2, 0~2)"));
    }
#endif // HUR_DEBUG

    real e = 0;
    if (i == 0 && j == 1 && k == 2) {
        e = 1;
    } else if (i == 1 && j == 2 && k == 0) {
        e = 1;
    } else if (i == 2 && j == 0 && k == 1) {
        e = 1;
    } else if (i == 0 && j == 2 && k == 1) {
        e = -1;
    } else if (i == 1 && j == 0 && k == 2) {
        e = -1;
    } else if (i == 2 && j == 1 && k == 0) {
        e = -1;
    }

    return e;
}