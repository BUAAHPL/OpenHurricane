/*!
 * \file timeMarching.cpp
 * \brief Main subroutines for time marching.
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

#include "timeMarching.hpp"
#include "EulerFlow.hpp"
#include "solver.hpp"

namespace OpenHurricane {
    createClassNameStr(timeMarching, "timeMarching");
    createObjFty(timeMarching, controller);
} // namespace OpenHurricane

void OpenHurricane::timeMarching::calcFluxSpectRadius() {
    const real prl = flow_.Prl();

    const faceList &faces = mesh().faces();
    const auto &fzl = mesh().faceZones();

    const auto &cells = mesh().cells();
    const auto &fA = mesh().faceArea();
    //const auto& cC = mesh().cellCentre();
    const auto &adjoinCCtr = mesh().adjoinCellCtr();
    const auto &cV = mesh().cellVolume();
    //const auto& fW = mesh().faceWeight();
    const auto &AR = mesh().aspectRatio();
    const auto &cellOrtho = mesh().cellOrtho();

    const auto &rho = flow_.rho();
    const auto &p = flow_.p();

    const auto &T = flow_.T();
    const auto &mu = flow_.mul();
    const auto &g = flow_.gama();

    for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
        if (fzl[fzi].isInterior() || fzl[fzi].isCutFace() || fzl[fzi].isPeriodic() ||
            fzl[fzi].isPeriodicShadow()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
                const auto &cl = faces[fi].leftCell();
                const auto &cr = faces[fi].rightCell();
                const real nn = fA[fi].magnitude();
                const vector &dd = adjoinCCtr[fi];
                const real rr = fabs(dd * fA[fi]) / nn;
                real ke = flow_.keEff(cl);
                real c = 0;
                real ffa = 1.0;
                if (solver_.useLowMachPrecon()) {
                    const real vv = v_[cl] * fA[fi].normalized();
                    ffa = 0.5 * (solver_.a4()[cl] + 1.0);
                    c = 0.5 *
                        sqrt(sqr(vv) * sqr(solver_.a4()[cl] - 1.0) + real(4) * solver_.a5()[cl]);
                } else {
                    c = sqrt(g[cl] * p[cl] / rho[cl]);
                }
                rai_[fi][0] = ffa * fabs(v_[cl] * fA[fi]) + c * nn;
                real mid;
                if (!isEuler_) {
                    mid = max(real(4.0 / 3.0), g[cl]) * ke / rho[cl];
                    rav_[fi][0] = mid * nn / rr;
                }

                ke = flow_.keEff(cr);
                c = 0;
                if (solver_.useLowMachPrecon()) {
                    const real vv = v_[cr] * fA[fi].normalized();
                    ffa = 0.5 * (solver_.a4()[cr] + 1.0);
                    c = 0.5 *
                        sqrt(sqr(vv) * sqr(solver_.a4()[cr] - 1.0) + real(4) * solver_.a5()[cr]);
                } else {
                    c = sqrt(g[cr] * p[cr] / rho[cr]);
                }
                rai_[fi][1] = ffa * fabs(v_[cr] * fA[fi]) + c * nn;
                if (!isEuler_) {
                    mid = max(real(4.0 / 3.0), g[cr]) * ke / rho[cr];
                    rav_[fi][1] = mid * nn / rr;
                }
                if (std::isnan(rav_[fi][0]) || std::isnan(rav_[fi][1])) {
                    LFatal("Spectral radius of viscous flux is NaN");
                }
            }
        } else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
                const auto &cl = faces[fi].leftCell();
                const auto &cr = faces[fi].rightCell();
                const real nn = fA[fi].magnitude();
                const vector &dd = adjoinCCtr[fi];
                const real rr = fabs(dd * fA[fi]) / nn;
                real ke = flow_.keEff(cl);
                real c = 0;
                real ffa = 1.0;
                if (solver_.useLowMachPrecon()) {
                    const real vv = v_[cl] * fA[fi].normalized();
                    ffa = 0.5 * (solver_.a4()[cl] + 1.0);
                    c = 0.5 *
                        sqrt(sqr(vv) * sqr(solver_.a4()[cl] - 1.0) + real(4) * solver_.a5()[cl]);
                } else {
                    c = sqrt(g[cl] * p[cl] / rho[cl]);
                }
                rai_[fi][0] = ffa * fabs(v_[cl] * fA[fi]) + c * nn;
                real mid;
                if (!isEuler_) {
                    mid = max(real(4.0 / 3.0), g[cl]) * ke / rho[cl];
                    rav_[fi][0] = mid * nn / rr;
                }

                rai_[fi][1] = 0.0;
                rav_[fi][1] = 0.0;
            }
        }
    }
}

OpenHurricane::timeMarching::timeMarching(const controller &cont, const runtimeMesh &mesh,
                                          const flowModel &flowM, solver &_solver,
                                          cellVectorArray &v)
    : mesh_(mesh), nPrimitives_(mesh.primitiveParamSize()), solver_(_solver), v_(v), flow_(flowM),
      objectList_(), countParam_(0), paramMap_(), dq_(mesh.nTotalCells()), dt_(_solver.dt()),
      dtReduceFactor_(0.5), shockFactor_(_solver.shockFactor()), isEuler_(false),
      rai_(object("rai", mesh, object::NOT_WRITE), mesh),
      rav_(object("rav", mesh, object::NOT_WRITE), mesh, Zero), pseudoTimePtr_(nullptr),
      JacPtr_(nullptr), sourceTermTreat_(sourceTermTreating::EXPLICITSOURCE) {
    if (cont.parent().found("flow")) {
        if (cont.parent().subController("flow").found("flowModel")) {
            string flowM = cont.parent().subController("flow").findWord("flowModel");
            if (flowM == EulerFlow::className_) {
                isEuler_ = true;
            }
        }
    }
    if (cont.found("timeMethod")) {
        const auto &timeCont = cont.subController("timeMethod");
        if (timeCont.found("sourceTermTreating")) {
            string sorw = timeCont.findWord("sourceTermTreating");
            trim(sorw);
            stringToUpperCase(sorw);
            if (sorw == "EXPLICITSOURCE") {
                sourceTermTreat_ = sourceTermTreating::EXPLICITSOURCE;
#ifdef HUR_DEBUG
                Pout << "    Info: using explicit source term treatment" << std::endl;
#endif // HUR_DEBUG
            } else if (sorw == "DIAGIMPLICISOURCE") {
                sourceTermTreat_ = sourceTermTreating::KIMDIAGIMPLICITSOURCE;
#ifdef HUR_DEBUG
                Pout << "    Info: using Kim diagonal Jacobian source term "
                        "treatment"
                     << std::endl;
#endif // HUR_DEBUG
            } else if (sorw == "FULLJACOBIAN") {
                sourceTermTreat_ = sourceTermTreating::FULLJACOBIAN;
#ifdef HUR_DEBUG
                Pout << "    Info: using full point-implicit Jacobian source "
                        "term treatment"
                     << std::endl;
#endif // HUR_DEBUG
            } else if (sorw == "FULLJACOBIANTABLE") {
                sourceTermTreat_ = sourceTermTreating::FULLJACOBIANTABLE;
#ifdef HUR_DEBUG
                Pout << "    Info: using full point-implicit Jacobian source "
                        "term treatment with tabulation"
                     << std::endl;
#endif // HUR_DEBUG
            } else {
                LFatal("Unknown string %s in %s", timeCont.findWord("sourceTermTreating").c_str(),
                       timeCont.name().c_str());
            }
        }
    }

    pseudoTimePtr_ =
        pseudoTime::creator(cont.subController("pseudoTime"), mesh, _solver.iter(), flow_,
                            _solver.dt(), _solver.shockFactor(), rai_, rav_,
                            _solver.temperatureFlag(), _solver.pressureFlag(), _solver.CFLFlag());
}

OpenHurricane::uniquePtr<OpenHurricane::timeMarching>
OpenHurricane::timeMarching::creator(const controller &cont, const runtimeMesh &mesh,
                                     const flowModel &flowM, solver &_solver, cellVectorArray &v) {
    if (!cont.found("timeMethod")) {
        LFatal("Cannot find timeMethod in %s", cont.name().c_str());
    }
    const auto &tmcont = cont.subController("timeMethod");

    string spschemeType = tmcont.findWord(timeMarching::className_);

    Pout << "    Info: setting time marching scheme: " << spschemeType << std::endl;
    defineInObjCreator(timeMarching, static_cast<std::string>(spschemeType), controller,
                       (cont, mesh, flowM, _solver, v));
}

OpenHurricane::timeMarching::~timeMarching() noexcept {
    for (integer i = 0; i < objectList_.size(); ++i) {
        objectList_[i] = nullptr;
    }
    pseudoTimePtr_.clear();
    JacPtr_.clear();
}

hur_nodiscard const OpenHurricane::iteration &OpenHurricane::timeMarching::iter() const noexcept {
    return solver_.iter();
}

hur_nodiscard OpenHurricane::real
OpenHurricane::timeMarching::vistParam(const integer paramI, const integer cellI) const {
    const integer objectI = paramMap_[paramI][0];
    object *ob = objectList_[objectI];
    if (ob->nElements() == 1) {
        const cellRealArray *f = static_cast<const cellRealArray *>(ob);

        return (*f)[cellI];
    } else if (ob->nElements() == 3) {
        const cellVectorArray *f = static_cast<const cellVectorArray *>(ob);

        return (*f)[cellI][paramMap_[paramI][1]];
    } else {
        return 0.0;
    }
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::timeMarching::vistParam(const integer cellI) const {
    realArray re(countParam_, Zero);

    for (integer pi = 0; pi < paramMap_.size(); ++pi) {
        re[pi] = vistParam(pi, cellI);
    }
    return re;
}

void OpenHurricane::timeMarching::pushParam(const integer paramI, const integer cellI,
                                            const real value) const {
    const integer objectI = paramMap_[paramI][0];
    object *ob = objectList_[objectI];
    if (ob->nElements() == 1) {
        cellRealArray *f = static_cast<cellRealArray *>(ob);

        (*f)[cellI] = value;
    } else if (ob->nElements() == 3) {
        cellVectorArray *f = static_cast<cellVectorArray *>(ob);

        (*f)[cellI][paramMap_[paramI][1]] = value;
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::timeMarching::vistRhs(const integer paramI,
                                                                       const integer cellI) const {
    const integer objectI = paramMap_[paramI][0];
    object *ob = objectList_[objectI];
    if (ob->nElements() == 1) {
        const cellRealArray *f = static_cast<const cellRealArray *>(ob);

        return (*f).rhs()[cellI];
    } else if (ob->nElements() == 3) {
        const cellVectorArray *f = static_cast<const cellVectorArray *>(ob);

        return (*f).rhs()[cellI][paramMap_[paramI][1]];
    } else {
        return 0.0;
    }
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::timeMarching::vistRhs(const integer cellI) const {
    realArray re(paramMap_.size(), Zero);

    for (integer pi = 0; pi < paramMap_.size(); ++pi) {
        re[pi] = vistRhs(pi, cellI);
    }
    return re;
}

hur_nodiscard OpenHurricane::real
OpenHurricane::timeMarching::vistDiagSource(const integer paramI, const integer cellI) const {
    const integer objectI = paramMap_[paramI][0];
    object *ob = objectList_[objectI];
    if (ob->nElements() == 1) {
        const cellRealArray *f = static_cast<const cellRealArray *>(ob);

        return (*f).diagSource()[cellI];
    } else if (ob->nElements() == 3) {
        const cellVectorArray *f = static_cast<const cellVectorArray *>(ob);

        return (*f).diagSource()[cellI][paramMap_[paramI][1]];
    } else {
        return 0.0;
    }
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::timeMarching::vistDiagSource(const integer cellI) const {
    realArray re(countParam_, Zero);

    for (integer pi = 0; pi < paramMap_.size(); ++pi) {
        re[pi] = vistDiagSource(pi, cellI);
    }
    return re;
}

void OpenHurricane::timeMarching::initializing() {
    for (integer ci = 0; ci < dq_.size(); ++ci) {
        dq_[ci].resize(countParam_, Zero);
    }
    if (fullJacobianSource()) {
        JacPtr_.reset(new cellRealSquareMatrixArray(mesh_));
        for (integer i = 0; i < JacPtr_->size(); ++i) {
            (*JacPtr_)[i].resize(countParam_);
            (*JacPtr_)[i].setZero();
        }
    } else if (fullJacobianSourceTable()) {
        JacPtr_.reset(new cellRealSquareMatrixArray(mesh_));
    }
}

void OpenHurricane::timeMarching::storeOldPrimitiveValue() {
    for (integer oi = 0; oi < objectList_.size(); ++oi) {
        object *ob = objectList_[oi];

        if (ob->nElements() == 1) {
            cellRealArray *f = static_cast<cellRealArray *>(ob);
            f->storeLastArray();
        } else if (ob->nElements() == 3) {
            cellVectorArray *f = static_cast<cellVectorArray *>(ob);
            f->storeLastArray();
        } else {
            errorAbortStr(("Number of components " + toString(ob->nElements()) +
                           " of parameter not supported yet"));
        }
    }
}

void OpenHurricane::timeMarching::timeStep() {
    calcFluxSpectRadius();
    calcShockFactor();
    if (iter().isSteadyFlow()) {
        if (solver_.hasTPLmtCells()) {
            pseudoTimes().cfl().setNotChange();
        }
    }
    pseudoTimes().restrictTimeStep();
}

void OpenHurricane::timeMarching::calcShockFactor() {
    const auto &p = flow_.p();
    const auto &cells = mesh().cells();
    const faceList &faces = mesh().faces();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        real ra_sum = 0.0;
        real ff = 1.0;
        for (integer i = 0; i < cells[n].faceSize(); ++i) {
            const integer j = cells[n].facei(i);
            const auto &cl = faces[j].leftCell();
            const auto &cr = faces[j].rightCell();

            ff = min(ff, min(p[cl], p[cr]) / (fabs(p[cr] - p[cl]) + tiny));
        }
        shockFactor_[n] = ff;
    }
    realTransfer myTransfer(mesh(), shockFactor_, false, true);
    myTransfer.transferInit();
    myTransfer.transferring();
}

hur_nodiscard bool OpenHurricane::timeMarching::isBDFScheme() const noexcept {
    return false;
}

hur_nodiscard bool OpenHurricane::timeMarching::isESDIRKScheme() const noexcept {
    return false;
}

OpenHurricane::pseudoTime &OpenHurricane::timeMarching::pseudoTimes() {
    return *pseudoTimePtr_;
}
