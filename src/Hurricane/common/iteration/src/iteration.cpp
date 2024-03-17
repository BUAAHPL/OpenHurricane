/*!
 * \file iteration.cpp
 * \brief The subroutines and functions of CFD time advance iteration
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
#include "iteration.hpp"
#include "flowModel.hpp"
#include "geometryMesh.hpp"
#include "hdf5I.hpp"
#include "hdf5O.hpp"
#include "monitors.hpp"
#include "runtimeMesh.hpp"
#include "solutionWrite.hpp"
#include "tecplotWriter.hpp"
#include "writeFaceZone.hpp"

#include "NullRefObj.hpp"

hur_nodiscard const OpenHurricane::iteration &OpenHurricane::iteration::nullObject() {
    return NullRefObj::nullRef<iteration>();
}

OpenHurricane::iteration::iteration(const char *_c, const argParse &arg)
    : registerTable(*this), cont_(_c), name_(_c), myMonitorPtr_(nullptr), solWritePtr_(nullptr),
      writeFaceZonePtr_(nullptr), isResConvergence_(false), cStep_(0), maxStep_(10000000),
      totalStep_(0), writeOutputStep_(500), restart_(false), restartFromUnsteady_(false),
      restartFrom_(), pTStepPtr_(nullptr), subIterPtr_(nullptr), refValues_(),
      flowState_(flowStateType::steady), timestepType_(timestepType::localTimeStep),
      readLastFromRelay_(false), writeLastToRelay_(false) {
    readCont();
}

OpenHurricane::iteration::iteration(const char *_c, const argParse &arg, const std::string &contStr)
    : registerTable(*this), cont_(_c), name_(_c), myMonitorPtr_(nullptr), solWritePtr_(nullptr),
      writeFaceZonePtr_(nullptr), isResConvergence_(false), cStep_(0), maxStep_(10000000),
      totalStep_(0), writeOutputStep_(500), restart_(false), restartFromUnsteady_(false),
      restartFrom_(), pTStepPtr_(nullptr), subIterPtr_(nullptr), refValues_(),
      flowState_(flowStateType::steady), timestepType_(timestepType::localTimeStep),
      readLastFromRelay_(false), writeLastToRelay_(false) {
    readCont(contStr);
}

OpenHurricane::iteration::~iteration() noexcept {
    clear();
}

hur_nodiscard const OpenHurricane::physicalTimeStep &OpenHurricane::iteration::pTStep() const {
    if (!hasPhysicalTimeStep()) {
        LFatal("Attempt to access a null physical time step pointer, in iteration: %s",
               name().c_str());
    }
    return *pTStepPtr_;
}

hur_nodiscard OpenHurricane::physicalTimeStep &OpenHurricane::iteration::pTStep() {
    if (!hasPhysicalTimeStep()) {
        LFatal("Attempt to access a null physical time step pointer, in iteration: %s",
               name().c_str());
    }
    return *pTStepPtr_;
}

hur_nodiscard const OpenHurricane::subIteration &OpenHurricane::iteration::subIter() const {
    if (!hasSubIteration()) {
        LFatal("Attempt to access a null sub-iteration pointer, in iteration: %s", name().c_str());
    }
    return *subIterPtr_;
}

hur_nodiscard bool OpenHurricane::iteration::end() const noexcept {
    bool runningIter = cStep_ >= maxStep_;
    if (hasPhysicalTimeStep()) {
        bool runningPhy = pTStepPtr_->end();
        return (runningIter || runningPhy);
    }
    return runningIter;
}

void OpenHurricane::iteration::clear() noexcept {
    cStep_ = 0;
    maxStep_ = 0;
    totalStep_ = 0;
    pTStepPtr_.clear();
    subIterPtr_.clear();

    myMonitorPtr_.clear();
    solWritePtr_.clear();
    writeFaceZonePtr_.clear();
}

hur_nodiscard bool OpenHurricane::iteration::iterating() noexcept {
    bool running = !end() && !isResConvergence_;
    if (running) {
        operator++();
    }
    return running;
}

void OpenHurricane::iteration::initializing(const flowModel &flows) {
    setMyMonitorPtr(new monitors(*this, flows.mesh()));
    setSolWrite(flows);

    setWriteFaceZonePtr(new writeFaceZone(flows, cont()));
    if (cont().found("ref")) {
        const auto &refCont = cont().subController("ref");
        refValues_.reset(refCont);
    }
}

OpenHurricane::iteration &OpenHurricane::iteration::operator++() noexcept {
    cStep_++;
    totalStep_++;
    if (hasPhysicalTimeStep()) {
        pTStepPtr_->operator++();
    }
    if (hasSubIteration()) {
        subIterPtr_->reset();
    }
    return *this;
}

