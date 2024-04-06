/*!
 * \file TVDRK4.cpp
 * \brief Main subroutines for 4th order TVDRK schemes.
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
#include "TVDRK4.hpp"
#include "solver.hpp"
namespace OpenHurricane {
    createClassNameStr(TVDRK4, "TVDRK4");
    registerObjFty(timeMarching, TVDRK4, controller);
} // namespace OpenHurricane

OpenHurricane::TVDRK4::TVDRK4(const controller &cont, const runtimeMesh &mesh,
                              const flowModel &flowM, solver &_solver, cellVectorArray &v)
    : timeMarching(cont, mesh, flowM, _solver, v), alpha_(4), beta_(4), stageSize_(4), rhoOi_(-1) {
    real cfl0 =
        cont.subController("timeMethod").subController("TVDRK4").findOrDefault<real>("cfl", 1.0);

    const real minCFL = min(pseudoTimes().cfl().CFLMax(), cfl0);

    pseudoTimes().cfl().setCFLMax(minCFL);
    pseudoTimes().unsetIsStretchAc();

    alpha_[0] = 1.0 / 2.0;
    alpha_[1] = 1.0 / 2.0;
    alpha_[2] = 1.0;
    alpha_[3] = 1.0 / 6.0;

    beta_[0] = 1;
    beta_[1] = 2;
    beta_[2] = 2;
    beta_[3] = 1;

    if (!explicitSource()) {
        LFatal("The TVD Runge-Kutta only supports explicit source terms treatment");
    }
}

OpenHurricane::TVDRK4::~TVDRK4() noexcept {}

void OpenHurricane::TVDRK4::stageUpdate() {
    storeOldPrimitiveValue();
    initOldResiduals();
    restoreOldResiduals(0);
    const auto &cV = mesh().cellVolume();
    cellRealArray *rho = static_cast<cellRealArray *>(objectList_[rhoOi_]);
    for (integer stageI = 0; stageI < stageSize_; ++stageI) {
        for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
            real dtdV = dt_[cellI] / cV[cellI];

            if (stageI == stageSize_ - 1) {
                for (integer oi = 0; oi < objectList_.size(); ++oi) {
                    object *ob = objectList_[oi];
                    if (oi == rhoOi_) {
                        (*rho)[cellI] = (*rho).lastArray()[cellI] +
                                        alpha_[stageI] * dtdV * (*rho).rhs().lastArray()[cellI];
                        continue;
                    }
                    if (ob->nElements() == 1) {
                        cellRealArray *f = static_cast<cellRealArray *>(ob);
                        (*f)[cellI] = (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] +
                                      alpha_[stageI] * dtdV * (*f).rhs().lastArray()[cellI];
                    } else if (ob->nElements() == 3) {
                        cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                        (*f)[cellI] = (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] +
                                      alpha_[stageI] * dtdV * (*f).rhs().lastArray()[cellI];
                    }
                }
            } else {
                for (integer oi = 0; oi < objectList_.size(); ++oi) {
                    object *ob = objectList_[oi];
                    if (oi == rhoOi_) {
                        (*rho)[cellI] =
                            (*rho).lastArray()[cellI] + alpha_[stageI] * dtdV * (*rho).rhs()[cellI];
                        continue;
                    }
                    if (ob->nElements() == 1) {
                        cellRealArray *f = static_cast<cellRealArray *>(ob);
                        (*f)[cellI] = (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] +
                                      alpha_[stageI] * dtdV * (*f).rhs()[cellI];
                    } else if (ob->nElements() == 3) {
                        cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                        (*f)[cellI] = (*rho).lastArray()[cellI] * (*f).lastArray()[cellI] +
                                      alpha_[stageI] * dtdV * (*f).rhs()[cellI];
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
            restoreOldResiduals(stageI + 1);
        }
    }
    rho = nullptr;
}

void OpenHurricane::TVDRK4::updatePrimitives() {
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

void OpenHurricane::TVDRK4::restoreOldResiduals(const integer stagei) {
    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];

        if (ob->nElements() == 1) {
            cellRealArray *f = static_cast<cellRealArray *>(ob);
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                f->rhs().lastArray()[cellI] += beta_[stagei] * f->rhs()[cellI];
            }
        } else if (ob->nElements() == 3) {
            cellVectorArray *f = static_cast<cellVectorArray *>(ob);
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                f->rhs().lastArray()[cellI] += beta_[stagei] * f->rhs()[cellI];
            }
        }
    }
}

void OpenHurricane::TVDRK4::initOldResiduals() {
    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];

        if (ob->nElements() == 1) {
            cellRealArray *f = static_cast<cellRealArray *>(ob);
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                f->rhs().lastArray()[cellI] = 0;
            }
        } else if (ob->nElements() == 3) {
            cellVectorArray *f = static_cast<cellVectorArray *>(ob);
            for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
                f->rhs().lastArray()[cellI] = 0;
            }
        }
    }
}

void OpenHurricane::TVDRK4::timeStep() {
    timeMarching::timeStep();
}

void OpenHurricane::TVDRK4::initializing() {
    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];
        if (ob->name() == "rho") {
            rhoOi_ = oi;
            break;
        }
    }
}

void OpenHurricane::TVDRK4::marching() {
    stageUpdate();
}
