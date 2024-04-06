/*!
 * \file TVDRK.cpp
 * \brief Main subroutines for TVDRK schemes.
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
#include "TVDRK.hpp"
#include "solver.hpp"
namespace OpenHurricane {
    createClassNameStr(TVDRK, "TVDRK");
    registerObjFty(timeMarching, TVDRK, controller);
} // namespace OpenHurricane

OpenHurricane::TVDRK::TVDRK(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
                            solver &_solver, cellVectorArray &v)
    : timeMarching(cont, mesh, flowM, _solver, v), stageSize_(2) {
    string stageCoefW = cont.subController("timeMethod")
                            .subController("TVDRK")
                            .findWordOrDefault("TVDRKOrder", "TVDRK2");
    setStageCoeffs(stageCoefW);
    real cfl0 =
        cont.subController("timeMethod").subController("TVDRK").findOrDefault<real>("cfl", 1.0);

    const real minCFL = min(pseudoTimes().cfl().CFLMax(), cfl0);

    pseudoTimes().cfl().setCFLMax(minCFL);
    pseudoTimes().unsetIsStretchAc();

    if (!explicitSource()) {
        LFatal("The TVD Runge-Kutta only supports explicit source terms treatment");
    }
}

OpenHurricane::TVDRK::~TVDRK() noexcept {}

void OpenHurricane::TVDRK::stageUpdate() {
    storeOldPrimitiveValue();
    auto &cV = mesh().cellVolume();
    cellRealArray *rho = static_cast<cellRealArray *>(objectList_[rhoOi_]);
    for (integer stageI = 0; stageI < stageSize_; ++stageI) {
        for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
            real dtdV = dt_[cellI] / cV[cellI];

            if (stageI == 0) {
                for (integer oi = 0; oi < objectList_.size(); ++oi) {
                    object *ob = objectList_[oi];
                    if (oi == rhoOi_) {
                        (*rho)[cellI] = (*rho).lastArray()[cellI] + dtdV * (*rho).rhs()[cellI];
                        continue;
                    }

                    if (ob->nElements() == 1) {
                        cellRealArray *f = static_cast<cellRealArray *>(ob);
                        (*f)[cellI] = (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] +
                                      dtdV * (*f).rhs()[cellI];
                    } else if (ob->nElements() == 3) {
                        cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                        (*f)[cellI] = (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] +
                                      dtdV * (*f).rhs()[cellI];
                    }
                }
            } else {
                for (integer oi = 0; oi < objectList_.size(); ++oi) {
                    object *ob = objectList_[oi];

                    if (oi == rhoOi_) {
                        (*rho)[cellI] = alpha_[stageI - 1][0] * (*rho).lastArray()[cellI] +
                                        alpha_[stageI - 1][1] * (*rho)[cellI] +
                                        alpha_[stageI - 1][2] * dtdV * (*rho).rhs()[cellI];
                        continue;
                    }

                    if (ob->nElements() == 1) {
                        cellRealArray *f = static_cast<cellRealArray *>(ob);
                        (*f)[cellI] = alpha_[stageI - 1][0] * (*f).lastArray()[cellI] *
                                          (*rho).lastArray()[cellI] +
                                      alpha_[stageI - 1][1] * (*f)[cellI] +
                                      alpha_[stageI - 1][2] * dtdV * (*f).rhs()[cellI];
                    } else if (ob->nElements() == 3) {
                        cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                        (*f)[cellI] = alpha_[stageI - 1][0] * (*rho).lastArray()[cellI] *
                                          (*f).lastArray()[cellI] +
                                      alpha_[stageI - 1][1] * (*f)[cellI] +
                                      alpha_[stageI - 1][2] * dtdV * (*f).rhs()[cellI];
                    }
                }
            }
        }

        updatePrimitives();
        if (stageI < stageSize_ - 1) {
            solver_.updatePrimitives(true);

            solver_.calculateFc();
            solver_.calculateSource();
            solver_.calculateFv();
            restoreConserves();
        }
    }
    rho = nullptr;
}

void OpenHurricane::TVDRK::updatePrimitives() {
    cellRealArray *rho = static_cast<cellRealArray *>(objectList_[rhoOi_]);
    for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
        for (integer oi = 0; oi < objectList_.size(); ++oi) {
            object *ob = objectList_[oi];
            if (oi == rhoOi_) {
                continue;
            }

            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f)[cellI] = (*f)[cellI] / (*rho)[cellI];
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f)[cellI] = (*f)[cellI] / (*rho)[cellI];
            }
        }
    }
}

void OpenHurricane::TVDRK::restoreConserves() {
    cellRealArray *rho = static_cast<cellRealArray *>(objectList_[rhoOi_]);
    for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
        for (integer oi = 0; oi < objectList_.size(); ++oi) {
            object *ob = objectList_[oi];
            if (oi == rhoOi_) {
                continue;
            }

            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                (*f)[cellI] = (*f)[cellI] * (*rho)[cellI];
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                (*f)[cellI] = (*f)[cellI] * (*rho)[cellI];
            }
        }
    }
}

void OpenHurricane::TVDRK::setStageCoeffs(const string &w) {
    if (w == "TVDRK2") {
        alpha_.resize(1);
        alpha_[0].resize(3);
        alpha_[0][0] = 0.5;
        alpha_[0][1] = 0.5;
        alpha_[0][2] = 0.5;
        stageSize_ = 2;
    } else if (w == "TVDRK3") {
        alpha_.resize(2);
        alpha_[0].resize(3);
        alpha_[1].resize(3);
        alpha_[0][0] = 0.75;
        alpha_[0][1] = 0.25;
        alpha_[0][2] = 0.25;
        alpha_[1][0] = 1.0 / 3.0;
        alpha_[1][1] = 2.0 / 3.0;
        alpha_[1][2] = 2.0 / 3.0;
        stageSize_ = 3;
    }
}

void OpenHurricane::TVDRK::timeStep() {
    timeMarching::timeStep();
}

void OpenHurricane::TVDRK::initializing() {
    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];
        if (ob->name() == "rho") {
            rhoOi_ = oi;
            break;
        }
    }
}

void OpenHurricane::TVDRK::marching() {
    stageUpdate();
}
