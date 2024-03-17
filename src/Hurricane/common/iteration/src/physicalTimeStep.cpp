/*!
 * \file physicalTimeStep.cpp
 * \brief The subroutines and functions of physical time.
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
#include "physicalTimeStep.hpp"

void OpenHurricane::physicalTimeStep::updateLastTimeStep() noexcept {
    for (integer i = lastTimeStep_.size() - 1; i > 0; --i) {
        lastTimeStep_[i] = lastTimeStep_[i - 1];
    }
    lastTimeStep_[0] = pTimeStep_;
}

OpenHurricane::physicalTimeStep::physicalTimeStep(const real timeStep, const real maxTime)
    : pTimeStep_(timeStep), maxTime_(maxTime), totalTime_(Zero), isDynamicTimeStep_(false),
      isDynamicCFLSet_(false), dyCFL_(0.8), lastTimeStep_() {
    lastTimeStep_.resize(3, Zero);
}

OpenHurricane::physicalTimeStep::physicalTimeStep(const real timeStep, const real maxTime,
                                              const real totalTime)
    : pTimeStep_(timeStep), maxTime_(maxTime), totalTime_(totalTime), isDynamicTimeStep_(false),
      isDynamicCFLSet_(false), dyCFL_(0.8), lastTimeStep_() {
    lastTimeStep_.resize(3, Zero);
}

void OpenHurricane::physicalTimeStep::setTimeStep(const real tS) noexcept {
    if (totalTime_ < maxTime_ && totalTime_ + tS > maxTime_) {
        pTimeStep_ = maxTime_ - totalTime_;
    } else {
        pTimeStep_ = tS;
    }
}

hur_nodiscard bool OpenHurricane::physicalTimeStep::timeMarching() noexcept {
    bool isEnd = end();
    if (!isEnd) {
        operator++();
    }
    return isEnd;
}

OpenHurricane::physicalTimeStep &OpenHurricane::physicalTimeStep::operator++() noexcept {
    if (totalTime_ < maxTime_ && (totalTime_ + pTimeStep_ > maxTime_)) {
        pTimeStep_ = maxTime_ - totalTime_;
        totalTime_ = maxTime_;
        updateLastTimeStep();
        return *this;
    }
    totalTime_ = totalTime_ + pTimeStep_;
    updateLastTimeStep();
    return *this;
}
