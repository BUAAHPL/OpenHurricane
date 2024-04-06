/*!
 * \file MRK.cpp
 * \brief Main subroutines for MRK schemes.
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
#include "MRK.hpp"
#include "solver.hpp"
namespace OpenHurricane {
    createClassNameStr(MRK, "MRK");
    registerObjFty(timeMarching, MRK, controller);
} // namespace OpenHurricane

OpenHurricane::MRK::MRK(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
                        solver &_solver, cellVectorArray &v)
    : timeMarching(cont, mesh, flowM, _solver, v), alpha_(), sigma_(), rhoOi_(-1),
      ep_(cont.subController("timeMethod").subController("MRK").findOrDefault<real>("ep", 0.5)),
      maxJacStep_(cont.subController("timeMethod")
                      .subController("MRK")
                      .findOrDefault<integer>("maxJacStep", 2)),
      cflRatioForImpResSmoo_(cont.subController("timeMethod")
                                 .subController("MRK")
                                 .findOrDefault<real>("cflRatioForImpResSmoo", 2.0)),
      impResSmooth_(false) {
    string stageCoefW =
        cont.subController("timeMethod").subController("MRK").findWordOrDefault("MRKStage", "f3");

    setStageCoeffs(stageCoefW);

    const auto &MRKCont = cont.subController("timeMethod").subController("MRK");
    impResSmooth_ = controllerSwitch(MRKCont)("impResSmooth", impResSmooth_);

    pseudoTimes().unsetIsStretchAc();

    if (!explicitSource()) {
        LFatal("The Multi-stage Runge-Kutta only supports explicit "
               "source terms treatment");
    }
}

OpenHurricane::MRK::~MRK() noexcept {}

void OpenHurricane::MRK::stageUpdate() {
    storeOldPrimitiveValue();
    auto &cV = mesh().cellVolume();
    cellRealArray *rho = static_cast<cellRealArray *>(objectList_[rhoOi_]);
    for (integer stageI = 0; stageI < alpha_.size(); ++stageI) {
        implicitResidualSmoothing();
        for (integer cellI = 0; cellI < mesh().nCells(); ++cellI) {
            real dtdV = dt_[cellI] / cV[cellI];

            for (integer oi = 0; oi < objectList_.size(); ++oi) {
                if (oi == rhoOi_) {
                    (*rho)[cellI] =
                        (*rho).lastArray()[cellI] + alpha_[stageI] * dtdV * (*rho).rhs()[cellI];
                    continue;
                }
                object *ob = objectList_[oi];

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

        updatePrimitives();
        if (stageI < alpha_.size() - 1) {
            solver_.updatePrimitives(true);
            solver_.calculateFc();
            solver_.calculateSource();
            solver_.calculateFv();
        }
    }
}

void OpenHurricane::MRK::updatePrimitives() {
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

void OpenHurricane::MRK::restoreConserves() {
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

void OpenHurricane::MRK::setStageCoeffs(const string &w) {
    if (w == "f3") {
        alpha_.resize(3);
        alpha_[0] = 0.1481;
        alpha_[1] = 0.400;
        alpha_[2] = 1.0;
        sigma_ = 1.5;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else if (w == "f4") {
        alpha_.resize(4);
        alpha_[0] = 0.0833;
        alpha_[1] = 0.2069;
        alpha_[2] = 0.4265;
        alpha_[3] = 1.0;
        sigma_ = 2.0;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else if (w == "f5") {
        alpha_.resize(5);
        alpha_[0] = 0.0533;
        alpha_[1] = 0.1263;
        alpha_[2] = 0.2375;
        alpha_[3] = 0.4414;
        alpha_[4] = 1.0;
        sigma_ = 2.5;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else if (w == "s3") {
        alpha_.resize(3);
        alpha_[0] = 0.1918;
        alpha_[1] = 0.4929;
        alpha_[2] = 1.0;
        sigma_ = 0.69;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else if (w == "s4") {
        alpha_.resize(4);
        alpha_[0] = 0.1084;
        alpha_[1] = 0.2602;
        alpha_[2] = 0.5052;
        alpha_[3] = 1.0;
        sigma_ = 0.92;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else if (w == "s5") {
        alpha_.resize(5);
        alpha_[0] = 0.0695;
        alpha_[1] = 0.1602;
        alpha_[2] = 0.2898;
        alpha_[3] = 0.5060;
        alpha_[4] = 1.0;
        sigma_ = 1.15;
        //cfl_ = sigma_;
        //pseudoTime_.cfl().setCFLMax(sigma_);
    } else {
        alpha_.resize(3);
        alpha_[0] = 0.1918;
        alpha_[1] = 0.4929;
        alpha_[2] = 1.0;
        sigma_ = 0.69;
    }
    pseudoTimes().cfl().setCFLMax(min(sigma_, pseudoTimes().cfl().CFLMax()));
}

void OpenHurricane::MRK::timeStep() {
    timeMarching::timeStep();
    if (impResSmooth_) {
        for (integer n = 0; n < mesh().nCells(); ++n) {
            dt_[n] *= cflRatioForImpResSmoo_;
        }
    }
}

void OpenHurricane::MRK::initializing() {
    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];
        if (ob->name() == "rho") {
            rhoOi_ = oi;
            break;
        }
    }
}

void OpenHurricane::MRK::marching() {
    stageUpdate();
}

void OpenHurricane::MRK::implicitResidualSmoothing(cellRealArray &Q) const {
    realArray Rm(Q.mesh().nTotalCells(), Zero); // R_m
    // R_(m-1)
    cellRealArray Rmm1(object(Q.name() + "_Rmm1", Q.mesh(), object::NOT_WRITE, object::TEMPORARY),
                       Q.mesh(), Zero);

    for (integer i = 0; i < mesh().nCells(); ++i) {
        Rmm1[i] = Q.rhs()[i];
    }

    const auto &faces = mesh().faces();
    const auto &cells = mesh().cells();
    for (integer m = 0; m < maxJacStep_; ++m) {
        //fv::transfer(Rmm1, true);
        realTransfer myTransfer(mesh(), Rmm1, false, true);
        myTransfer.transferInit();
        myTransfer.transferring();
        for (integer i = 0; i < mesh().nCells(); ++i) {
            const auto &fl = cells[i].facesList();

            real Rjm = Zero;
            integer Nf = fl.size();
            for (integer fli = 0; fli < fl.size(); ++fli) {
                const integer fi = fl[fli];
                const auto &cl = faces[fi].leftCell();
                const auto &cr = faces[fi].rightCell();

                const integer j = cl + cr - i;

                if (j < mesh().nCells() || (j >= mesh().nCells() && mag(Rmm1[j]) > tiny)) {
                    Nf--;
                    Rjm += Rmm1[j];
                }
            }

            Rm[i] = (Q.rhs()[i] + ep_ * Rjm) / (real(1.0) + ep_ * real(Nf));
        }

        Rmm1 = Rm;
    }

    for (integer i = 0; i < Rm.size(); ++i) {
        Q.rhs()[i] = Rm[i];
    }
}

void OpenHurricane::MRK::implicitResidualSmoothing(cellVectorArray &Q) const {
    vectorArray Rm(mesh().nCells()); // R_m

    // R_(m-1)
    cellVectorArray Rmm1(object(Q.name() + "_Rmm1", Q.mesh(), object::NOT_WRITE, object::TEMPORARY),
                         Q.mesh(), Zero);
    for (integer i = 0; i < mesh().nCells(); ++i) {
        Rmm1[i] = Q.rhs()[i];
    }

    const auto &faces = mesh().faces();
    const auto &cells = mesh().cells();
    for (integer m = 0; m < maxJacStep_; ++m) {
        //fv::transfer(Rmm1, true);
        vectorTransfer myTransfer(mesh(), Rmm1, false, true);
        myTransfer.transferInit();
        myTransfer.transferring();
        for (integer i = 0; i < mesh().nCells(); ++i) {
            const auto &fl = cells[i].facesList();

            vector Rjm = Zero;
            integer Nf = fl.size();
            for (integer fli = 0; fli < fl.size(); ++fli) {
                const integer fi = fl[fli];
                const auto &cl = faces[fi].leftCell();
                const auto &cr = faces[fi].rightCell();

                const integer j = cl + cr - i;

                if (j < mesh().nCells() || (j >= mesh().nCells() && mag(Rmm1[j]) > tiny)) {
                    Nf--;
                    Rjm += Rmm1[j];
                }
            }

            Rm[i] = (Q.rhs()[i] + ep_ * Rjm) / (real(1.0) + ep_ * real(Nf));
        }

        Rmm1 = Rm;
    }

    for (integer i = 0; i < Rm.size(); ++i) {
        Q.rhs()[i] = Rm[i];
    }
}

void OpenHurricane::MRK::implicitResidualSmoothing() {
    if (!impResSmooth_) {
        return;
    }

    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];

        if (ob->nElements() == 1) {
            cellRealArray *Qptr = static_cast<cellRealArray *>(ob);
            implicitResidualSmoothing(*Qptr);
        } else if (ob->nElements() == 3) {
            cellVectorArray *Qptr = static_cast<cellVectorArray *>(ob);
            implicitResidualSmoothing(*Qptr);
        }
    }
}
