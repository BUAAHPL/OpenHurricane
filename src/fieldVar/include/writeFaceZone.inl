/*!
 * \file writeFaceZone.inl
 * \brief In-Line subroutines of the <i>writeFaceZone.hpp</i> file.
 * \author Chen Zhenyi
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

inline const OpenHurricane::stringList &
OpenHurricane::writeFaceZone::writeFaceZoneNameList() const noexcept {
    return writeFaceZoneNameList_;
}

inline const OpenHurricane::List<OpenHurricane::stringList> &
OpenHurricane::writeFaceZone::writeFaceZoneVarList() const noexcept {
    return writeFaceZoneVarList_;
}

inline void OpenHurricane::writeFaceZone::write() const {
    const auto &iter = flows_.mesh().Iteration();
    if (iter.hasPhysicalTimeStep() && firstCalling_) {
        const auto cTT = iter.pTStep().totalTime();
        timeElasped_ = 0;
        lastTime_ = cTT - iter.pTStep().pTimeStep();
        firstCalling_ = false;
    }
    if (writeNow() || samplingNow()) {
        updating();
    }

    if (writeNow()) {
        writeToFile();
    }
}

inline bool OpenHurricane::writeFaceZone::samplingNow() const {
    const auto &iter = flows_.mesh().Iteration();
    if (isSampling_ && iter.hasPhysicalTimeStep()) {
        if (iter.cStep() % samplingStep_ == 0) {
            return true;
        }
    }
    return false;
}

inline bool OpenHurricane::writeFaceZone::writeNow() const {
    if (writeFaceZoneNameList_.size() == 0) {
        return false;
    } else {
        const auto &iter = flows_.mesh().Iteration();
        if (iter.cStep() % iter.writeOutputStep() == 0 ||
            (iter.end() && !argParse::noWriteAtEnd())) {
            return true;
        }
    }
    return false;
}

hur_nodiscard inline OpenHurricane::faceZoneFieldParameterFuncMap &
OpenHurricane::writeFaceZone::faceZoneVarFuncMap() noexcept {
    return faceZoneVarFuncMap_;
}

hur_nodiscard inline const OpenHurricane::faceZoneFieldParameterFuncMap &
OpenHurricane::writeFaceZone::faceZoneVarFuncMap() const noexcept {
    return faceZoneVarFuncMap_;
}
