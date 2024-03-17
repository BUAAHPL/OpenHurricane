#include "iteration.hpp"
/*!
 * \file iteration.inl
 * \brief In-Line subroutines of the <i>iteration.hpp</i> file.
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

inline void OpenHurricane::iteration::setResConveg() const {
    isResConvergence_ = true;
}

hur_nodiscard inline bool OpenHurricane::iteration::restart() const noexcept {
    return restart_;
}

hur_nodiscard inline bool OpenHurricane::iteration::restartFromUnsteady() const noexcept {
    return restartFromUnsteady_;
}

hur_nodiscard inline const OpenHurricane::string &OpenHurricane::iteration::name() const noexcept {
    return name_;
}

hur_nodiscard inline OpenHurricane::fileName OpenHurricane::iteration::meshFileName() const {
    fileName outN = caseFileName();
    auto pathOut = outN.parentPath();
    auto fname = outN.name(true);
    string fext;
    fext = ".h5";
    fname += "Grid";
    outN = fname + fext;
    outN = pathOut / outN;
    return outN;
}

inline void OpenHurricane::iteration::changeName(const string &nN) {
    name_ = nN;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::iteration::cStep() const noexcept {
    return cStep_;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::iteration::maxStep() const noexcept {
    return maxStep_;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::iteration::totalStep() const noexcept {
    return totalStep_;
}

inline OpenHurricane::integer OpenHurricane::iteration::setTotalStep(const integer newTS) noexcept {
    integer tmpTS = totalStep_;
    totalStep_ = newTS;
    return tmpTS;
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::iteration::writeOutputStep() const noexcept {
    return writeOutputStep_;
}

hur_nodiscard inline bool OpenHurricane::iteration::hasPhysicalTimeStep() const noexcept {
    return !pTStepPtr_.isNull();
}

hur_nodiscard inline bool OpenHurricane::iteration::hasSubIteration() const noexcept {
    return !subIterPtr_.isNull();
}

hur_nodiscard inline bool OpenHurricane::iteration::subLoop()noexcept {
    if (hasSubIteration()) {
        return subIterPtr_->subIterating();
    } else {
        return false;
    }
}

inline void OpenHurricane::iteration::resetSubiter() {
    if (hasSubIteration()) {
        subIterPtr_->reset();
    }
}

hur_nodiscard inline OpenHurricane::integer OpenHurricane::iteration::leftStep() const noexcept {
    return maxStep_ - cStep_;
}

inline void OpenHurricane::iteration::setMaxTime(const real mT) noexcept {
    if (hasPhysicalTimeStep()) {
        pTStepPtr_->setMaxTime(mT);
    }
}

inline void OpenHurricane::iteration::setTimeStep(const real tS) noexcept {
    if (hasPhysicalTimeStep()) {
        pTStepPtr_->setTimeStep(tS);
    }
}

inline void OpenHurricane::iteration::setWriteLastToRelay(const bool isw) noexcept {
    writeLastToRelay_ = isw;
}

inline void OpenHurricane::iteration::setReadLastFromRelay(const bool isw) noexcept {
    readLastFromRelay_ = isw;
}

hur_nodiscard inline bool OpenHurricane::iteration::isReadLastFromRelay() const noexcept {
    return readLastFromRelay_;
}

hur_nodiscard inline bool OpenHurricane::iteration::isWriteLastToRelay() const noexcept {
    return writeLastToRelay_;
}

hur_nodiscard inline const OpenHurricane::referenceValues &
OpenHurricane::iteration::refValues() const noexcept {
    return refValues_;
}

hur_nodiscard inline OpenHurricane::referenceValues &OpenHurricane::iteration::refValues() noexcept {
    return refValues_;
}

hur_nodiscard inline bool OpenHurricane::iteration::isSteadyFlow() const noexcept {
    return flowState_ == flowStateType::steady;
}

hur_nodiscard inline bool OpenHurricane::iteration::isUnsteadyFlow() const noexcept {
    return flowState_ == flowStateType::unsteady;
}

hur_nodiscard inline bool OpenHurricane::iteration::isBeginLimit() const {
    if (!isUnsteadyFlow()) {
        return true;
    } else {
        if (subIterPtr_) {
            if (subIterPtr_->cSubStep() < (subIterPtr_->maxSubStep() / 2)) {
                return true;
            }
        }
        return false;
    }
}