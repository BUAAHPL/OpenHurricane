/*!
 * \file pseudoTime.cpp
 * \brief Main subroutines for time marching.
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

#include "pseudoTime.hpp"
namespace OpenHurricane {
    createClassNameStr(pseudoTime, "pseudoTime");
    createObjFty(pseudoTime, controller);
} // namespace OpenHurricane

OpenHurricane::pseudoTime::pseudoTime(const controller &cont, const runtimeMesh &mesh,
                                      const iteration &iter, const flowModel &flows, realArray &dt,
                                      const cellRealArray &shockFactor,
                                      const faceVector2DArray &rai, const faceVector2DArray &rav,
                                      const integerArray &temperatureFlag,
                                      const integerArray &pressureFlag,
                                      const cellIntegerArray &CFLFlag)
    : mesh_(mesh), iter_(iter), flow_(flows), dt_(dt), shockFactor_(shockFactor), rai_(rai),
      rav_(rav), temperatureFlag_(temperatureFlag), pressureFlag_(pressureFlag), CFLFlag_(CFLFlag),
      dtReduceFactor_(0.5), cflPtr_(nullptr),
      CForTimeStep_(cont.findOrDefault<real>("CForTimeStep", 1.0)), isStretchAc_(true),
      cflRatioMax_(cont.findOrDefault<real>("cflRatioMax", 1000.0)),
      minStretchScale_(cont.findOrDefault<real>("minStretchScale", 3.0)),
      minCellOrthogonality_(cont.findOrDefault<real>("minCellOrthogonality", 0.7)),
      isShockReduce_(true), minShockReduceFct_(0.3) {
    if (cont.found("isStretchAc")) {
        isStretchAc_ = controllerSwitch(cont)("isStretchAc", isStretchAc_);
        if (cont.found("StretchAc")) {
            const auto &strCont = cont.subController("StretchAc");
            minStretchScale_ = strCont.findOrDefault<real>("minStretchScale", 3.0);
            cflRatioMax_ = strCont.findOrDefault<real>("cflRatioMax", 1000.0);
        }
    }
    if (cont.found("CFLSetting")) {
        const auto &cflCont = cont.subController("CFLSetting");
        cflPtr_ = CFL::creator(iter_, mesh_, cflCont);
    } else {
        cflPtr_ = CFL::creator(iter_, mesh_, cont);
    }
    if (cont.found("shockReduce")) {
        const auto &srCont = cont.subController("shockReduce");
        controllerSwitch mySRContS(srCont);
        isShockReduce_ = mySRContS("isShockReduce", isShockReduce_);
        minShockReduceFct_ = srCont.findOrDefault<real>("minShockReduceFct", minShockReduceFct_);
    }
}

OpenHurricane::uniquePtr<OpenHurricane::pseudoTime> OpenHurricane::pseudoTime::creator(
    const controller &cont, const runtimeMesh &mesh, const iteration &iter, const flowModel &flows,
    realArray &dt, const cellRealArray &shockFactor, const faceVector2DArray &rai,
    const faceVector2DArray &rav, const integerArray &temperatureFlag,
    const integerArray &pressureFlag, const cellIntegerArray &CFLFlag) {
    string timeStepType = cont.findWord("timeStepMethod");
    defineInObjCreator(pseudoTime, timeStepType, controller,
                       (cont, mesh, iter, flows, dt, shockFactor, rai, rav, temperatureFlag,
                        pressureFlag, CFLFlag));
}

OpenHurricane::pseudoTime::~pseudoTime() noexcept {
    cflPtr_.clear();
}

void OpenHurricane::pseudoTime::limitFlagReduce() {
    for (integer n = 0; n < mesh_.nCells(); ++n) {
        if (temperatureFlag_[n] != 1 || pressureFlag_[n] != 1 || CFLFlag_[n] != 1) {
            dt_[n] *= dtReduceFactor_;
        }
    }
}

void OpenHurricane::pseudoTime::nonorthogonalMeshReduce() {
    const auto &cellOrtho = mesh_.cellOrtho();
    for (integer n = 0; n < mesh_.nCells(); ++n) {
        if (cellOrtho[n] < minCellOrthogonality_) {
            dt_[n] *= (max(cellOrtho[n] * real(0.5), real(0.1)));
        }
    }
}

void OpenHurricane::pseudoTime::stretchedMeshAccelerate() {
    stretchedMeshAccelerate(dt_);
}

void OpenHurricane::pseudoTime::stretchedMeshAccelerate(realArray &dt) {
    const auto &AR = mesh_.aspectRatio();
    for (integer n = 0; n < mesh_.nCells(); ++n) {
        if (temperatureFlag_[n] != 1 || pressureFlag_[n] != 1 || CFLFlag_[n] != 1) {
            continue;
        } else {
            if (!cflPtr_->isInitialStage()) {
                //Note AR is the cell aera aspect ratio. It is used to accelerate the convergence
                //for stretched meshes
                if (AR[n] > minStretchScale_) {
                    dt[n] *= min(real(cflRatioMax_), AR[n] / minStretchScale_);
                }
            }
        }
    }
}

OpenHurricane::real OpenHurricane::pseudoTime::globalTimeStep() const {
    return globalTimeStep(dt_);
}

OpenHurricane::real OpenHurricane::pseudoTime::globalTimeStep(realArray &dt) const {
    real dtmin = large;
    for (integer n = 0; n < mesh_.nCells(); ++n) {
        dtmin = min(dtmin, dt[n]);
    }
    HurMPI::allReduce(dtmin, MPI_MIN);
    dt = dtmin;
    return dtmin;
}

OpenHurricane::real OpenHurricane::pseudoTime::restrictTimeStep() {
    real dtmin = tiny;
    computingTimeStep();
    shockReduce();

    if (isStretchAc_) {
        stretchedMeshAccelerate();
    }
    nonorthogonalMeshReduce();
    //limitFlagReduce();

    if (iter_.isGlobalTimeStep()) {
        dtmin = globalTimeStep();
    }
    return dtmin;
}

OpenHurricane::real OpenHurricane::pseudoTime::timeStep() {
    real dtmin = tiny;
    computingTimeStep();
    if (iter_.isGlobalTimeStep()) {
        dtmin = globalTimeStep();
    }
    return dtmin;
}

OpenHurricane::real OpenHurricane::pseudoTime::getGlobalTimeStep(const real cfl0,
                                                                 const bool noScale) {
    real dtmin = tiny;
    realArray dtt(mesh_.nCells(), Zero);
    computingTimeStep(dtt, cfl0);

    if (!noScale) {
        if (isStretchAc_) {
            stretchedMeshAccelerate(dtt);
        }
    }

    dtmin = globalTimeStep(dtt);

    return dtmin;
}
