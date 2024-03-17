#include "ODEsSolver.hpp"
/*!
 * \file ODEsSolver.inl
 * \brief In-Line subroutines of the <i>ODEsSolver.hpp</i> file.
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

#ifdef TEST_PROCESS_TIME
hur_nodiscard inline OpenHurricane::integer OpenHurricane::ODEsSolver::countIter() const noexcept {
    return countIter_;
}
#endif

hur_nodiscard inline OpenHurricane::integer OpenHurricane::ODEsSolver::nEqs() const noexcept {
    return nEqs_;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::ODEsSolver::maxEqs() const noexcept {
    return maxEqns_;
}

inline bool OpenHurricane::ODEsSolver::reset(const integer nEqsNew) {
    if (nEqsNew != nEqs_) {
        nEqs_ = nEqsNew;
        if (nEqs_ > maxEqns_) {
            errorAbortStr(("Cannot set the number (" + toString(nEqs_) +
                              ") of equations larger than " +
                              toString(maxEqns_)));
        }
        resetArray(ATOL_);
        resetArray(RTOL_);
        resetArray(yTemp_);
        resetArray(dydt0_);
        return true;
    } else {
        return false;
    }
}

template <typename Type>
inline void OpenHurricane::ODEsSolver::resetArray(List<Type> &ar, const integer n) {
    ar.resize(n);
}

template <typename Type> inline void OpenHurricane::ODEsSolver::resetArray(List<Type> &ar) {
    resetArray(ar, nEqs_);
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::ODEsSolver::maxSteps() const noexcept {
    return maxSteps_;
}

hur_nodiscard inline const OpenHurricane::realArray &OpenHurricane::ODEsSolver::ATOL() const {
    return ATOL_;
}

hur_nodiscard inline const OpenHurricane::realArray &OpenHurricane::ODEsSolver::RTOL() const {
    return RTOL_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::ODEsSolver::safeScale() const noexcept {
    return safeScale_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::ODEsSolver::minScale() const noexcept {
    return minScale_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::ODEsSolver::maxScale() const noexcept {
    return maxScale_;
}

hur_nodiscard inline OpenHurricane::real OpenHurricane::ODEsSolver::alphaInc() const noexcept {
    return alphaInc_;
}