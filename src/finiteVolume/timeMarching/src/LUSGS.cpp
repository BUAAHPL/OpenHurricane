/*!
 * \file LUSGS.cpp
 * \brief Main subroutines for LUSGS schemes.
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
#include "LUSGS.hpp"
#include "solver.hpp"
namespace OpenHurricane {
    createClassNameStr(LUSGS, "LUSGS");
    registerObjFty(timeMarching, LUSGS, controller);
} // namespace OpenHurricane

OpenHurricane::LUSGS::LUSGS(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
                            solver &_solver, cellVectorArray &v, const bool modifiedDiagnalTime)
    : timeMarching(cont, mesh, flowM, _solver, v),
      omega_(cont.subController("timeMethod")
                 .subController("LUSGS")
                 .findOrDefault<real>("omegaForLUSGS", 1.1)),
      betaT_(cont.subController("timeMethod")
                 .subController("LUSGS")
                 .findOrDefault<real>("betaT", 0.2)),
      modifiedDiagnalTime_(modifiedDiagnalTime), alphaDT_(0.0), scaleTimeStep_(1.0),
      fullJacFlagPtr_(nullptr) {}

OpenHurricane::LUSGS::~LUSGS() noexcept {
    HurDelete(fullJacFlagPtr_);
}

void OpenHurricane::LUSGS::lowerLoop() {
    const integer nTotalCells = mesh().nTotalCells();
    const integer nCells = mesh().nCells();
    const auto &cV = mesh().cellVolume();
    const auto &fA = mesh().faceArea();
    const auto &cC = mesh().cellCentre();
    const auto &faceL = mesh().faces();
    const auto &cellL = mesh().cells();
    auto &rho = const_cast<flowModel &>(flow_).rho();
    for (integer ci = 0; ci < nTotalCells; ++ci) {
        dq_[ci] = Zero;
    }

    if (fullJacobianSource()) {
        fullJacFlag() = false;
    }

    realArray flux(countParam());
    realArray adq(countParam());
    for (integer n = 0; n < nCells; ++n) {
        real dm;
        dm = cV[n] / (scaleTimeStep_ * dt_[n]);
        if (modifiedDiagnalTime_) {
            dm += alphaDT_ * cV[n];
        }
        flux = Zero;
        adq = Zero;

        for (integer fi = 0; fi < cellL[n].faceSize(); fi++) {
            const integer faceI = cellL[n].facei(fi);
            const auto &cl = faceL[faceI].leftCell();
            const auto &cr = faceL[faceI].rightCell();

            //const integer m = faceL[faceI].oppositeCell(n);
            const integer m = (n == cl) ? cr : cl;

            real ram;
            real bsign;

            if (n == cl) {
                dm += 0.5 * omega_ * rai_[faceI][0] + rav_[faceI][0];
                ram = 0.5 * rai_[faceI][1] + rav_[faceI][1];
            } else {
                dm += 0.5 * omega_ * rai_[faceI][1] + rav_[faceI][1];
                ram = 0.5 * rai_[faceI][0] + rav_[faceI][0];
            }
            if (m < n) {
                if (n == cl) {
                    bsign = -1.0;
                } else {
                    bsign = 1.0;
                }

                fluxJacoDQ(m, bsign * fA[faceI], adq);
                for (integer ii = 0; ii < flux.size(); ++ii) {
                    flux[ii] += (0.5 * adq[ii] - ram * dq_[m][ii]);
                }
            }
        }

        for (integer pi = 0; pi < paramMap_.size(); ++pi) {
            flux[pi] = vistRhs(pi, n) - flux[pi];
        }
        //flux = vistRhs(n) - flux;

        if (explicitSource()) {
            for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                dq_[n][pi] = flux[pi] / dm;
            }
        } else if (diagonalImpSource()) {
            for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                dq_[n][pi] = flux[pi] / (dm - vistDiagSource(pi, n));
            }
        } else if (fullJacobianSource()) {
            Jacobian()[n] *= -1;
            for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                Jacobian()[n](pi, pi) = dm + Jacobian()[n](pi, pi);
            }

            //Jacobian()[n].inverseAndStore();
            Jacobian()[n] = inv(Jacobian()[n]);

            dq_[n] = Jacobian()[n] * flux;

            if (std::isnan(dq_[n][0])) {
                if (!std::isnan(flux[0])) {
                    fullJacFlag()[n] = true;
                    real dmm = dm - cV[n] / (scaleTimeStep_ * dt_[n]);
                    dmm += cV[n] / (0.1 * scaleTimeStep_ * dt_[n]);
                    for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                        dq_[n][pi] = flux[pi] / dm;
                    }
                }
            }

        } else {
            LFatal("Unsupported source term treatment");
        }

        if (std::isnan(dq_[n][0])) {
            std::cout << "(Lower-looping) dq0 is not a number in cell: " << n
                      << ", cell centre = " << toString(cC[n]) << std::endl;
            std::cout << "rhs = " << vistRhs(n)[0] << ", dm = " << dm << ", dt = " << dt_[n]
                      << ", flux[0] = " << flux[0] << std::endl;

            HurMPI::abort(HurMPI::getComm(), EXIT_FAILURE);
        }

        if (std::isinf(dq_[n][1])) {
            std::cout << "(Lower-looping) dq1 is infinite in cell: " << n
                      << ", cell centre = " << toString(cC[n]) << std::endl;
            std::cout << "rhs = " << vistRhs(n)[1] << ", dm = " << dm << ", dt = " << dt_[n]
                      << ", flux[1] = " << flux[1] << ", adq[1] = " << adq[1] << std::endl;
            HurMPI::abort(HurMPI::getComm(), EXIT_FAILURE);
        }
        if (std::isnan(dq_[n][4])) {
            auto &T = flow_.T();
            auto &E = flow_.E();
            std::cout << "(Lower-looping) dq4 is not a number in cell: " << n << std::endl;
            std::cout << "rhs = " << vistRhs(n)[4] << ", T = " << T[n] << ", E = " << E[n]
                      << ", flux[4] = " << flux[4] << std::endl;
            HurMPI::abort(HurMPI::getComm(), EXIT_FAILURE);
        }
    }
}

void OpenHurricane::LUSGS::upperLoop() {
    const integer nTotalCells = mesh().nTotalCells();
    const integer nCells = mesh().nCells();
    const auto &cV = mesh().cellVolume();
    const auto &fA = mesh().faceArea();
    const auto &cC = mesh().cellCentre();
    const auto &faceL = mesh().faces();
    const auto &cellL = mesh().cells();
    realArray flux(countParam());
    realArray adq(countParam());

    for (integer n = nCells - 1; n >= 0; --n) {
        flux = Zero;
        adq = Zero;
        real dm = cV[n] / (scaleTimeStep_ * dt_[n]);
        if (modifiedDiagnalTime_) {
            dm += alphaDT_ * cV[n];
        }
        for (integer fi = 0; fi < cellL[n].faceSize(); fi++) {
            const integer faceI = cellL[n].facei(fi);
            const auto &cl = faceL[faceI].leftCell();
            const auto &cr = faceL[faceI].rightCell();

            const integer m = faceL[faceI].oppositeCell(n);
            //const integer m = (n == cl) ? cr : cl;

            real ram;
            real bsign;

            if (n == cl) {
                dm += 0.5 * omega_ * rai_[faceI][0] + rav_[faceI][0];
                ram = 0.5 * rai_[faceI][1] + rav_[faceI][1];
            } else {
                dm += 0.5 * omega_ * rai_[faceI][1] + rav_[faceI][1];
                ram = 0.5 * rai_[faceI][0] + rav_[faceI][0];
            }
            if (m > n) {
                if (n == cl) {
                    bsign = -1.0;
                } else {
                    bsign = 1.0;
                }
                fluxJacoDQ(m, bsign * fA[faceI], adq);
                for (integer ii = 0; ii < flux.size(); ++ii) {
                    flux[ii] += (0.5 * adq[ii] - ram * dq_[m][ii]);
                }
            }
        }

        if (explicitSource()) {
            for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                dq_[n][pi] -= flux[pi] / dm;
            }
        } else if (diagonalImpSource()) {
            for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                dq_[n][pi] -= flux[pi] / (dm - vistDiagSource(pi, n));
            }
        } else if (fullJacobianSource()) {
            if (fullJacFlag()[n]) {
                for (integer pi = 0; pi < paramMap_.size(); ++pi) {
                    real dmm = dm - cV[n] / (scaleTimeStep_ * dt_[n]);
                    dmm += cV[n] / (0.1 * scaleTimeStep_ * dt_[n]);
                    dq_[n][pi] -= flux[pi] / dmm;
                }
            } else {
                dq_[n] -= (Jacobian()[n] * flux);
            }
        } else {
            LFatal("Unsupported source term treatment");
        }
    }

    if (fullJacobianSource()) {
        for (integer n = 0; n < nCells; ++n) {
            Jacobian()[n].setZero();
        }
    }
}

void OpenHurricane::LUSGS::update() {
    const real coef1 = -0.20;
    const real coef2 = 2.0;

    /**\brief Last time-step value of primitive variables.*/
    realArray qold(countParam());
    auto &rho = const_cast<flowModel &>(flow_).rho();
    auto &p = const_cast<flowModel &>(flow_).p();
    auto &T = const_cast<flowModel &>(flow_).T();
    auto &E = const_cast<flowModel &>(flow_).E();
    auto &iflag = const_cast<flowModel &>(flow_).temperatureFlag();
    //auto& mixtures = const_cast<flowModel&>(flow_).mixtures();
    real sumRelCorrec = Zero;

    integer countNan = 0;

    for (integer n = 0; n < mesh().nCells(); ++n) {
        qold[0] = rho[n];
        qold[1] = v_[n][0];
        qold[2] = v_[n][1];
        qold[3] = v_[n][2];
        qold[4] = E[n];
        for (integer j = 5; j < countParam(); j++) {
            qold[j] = vistParam(j, n);
        }
        real sigmaT = 1.0;
        integer countMT = 0;

        while (true) {
            real rhon = qold[0];
            real rhoOld = rhon;
            real rhou = rhon * qold[1];
            real rhov = rhon * qold[2];
            real rhow = rhon * qold[3];
            real rhoE = rhon * qold[4];

            /**\brief New value of conservative variables.*/
            rhon += sigmaT * dq_[n][0];
            rhou += sigmaT * dq_[n][1];
            rhov += sigmaT * dq_[n][2];
            rhow += sigmaT * dq_[n][3];
            rhoE += sigmaT * dq_[n][4];
            real de = dq_[n][4];
            if (rhon < 0.0) {
                pseudoTimes().cfl().setBreakdowns();
                break;
            }
            rho[n] = rhon;
            v_[n][0] = rhou / rhon;
            v_[n][1] = rhov / rhon;
            v_[n][2] = rhow / rhon;
            E[n] = rhoE / rhon;

            for (integer j = 5; j < countParam(); j++) {
                real qOld = qold[j];
                real rhoq = rhoOld * qOld;
                rhoq += sigmaT * dq_[n][j];
                real qnew = rhoq / rhon;
                if (isnan(qnew)) {
                    qnew = qOld;
                    countNan++;
                }
                //qnew = max(min(qnew, real((1.0 + 0.05)) * qold[j]), real((1.0 - 0.05)) * qold[j]);
                pushParam(j, n, qnew);
            }

            // New temperature
            real newT =
                flow_.mixtures().thermalTable().TEa_rho((E[n] - real(0.5) * v_[n].magSqr()), rho[n],
                                                        T[n], iflag[n], flow_.mixtures().Yi(), n);
            real Td = T[n];
            real deltaT = mag(Td - newT);
            if (deltaT > betaT_ * Td && countMT < 2) {
                countMT++;
                sigmaT = (betaT_ * Td / deltaT) * sigmaT;
            } else {
                T[n] = newT;
                sumRelCorrec += sqr(sigmaT * dq_[n][0] / rhoOld);
                break;
            }
        }

        if (iflag[n] == 0 || iflag[n] == 2) {

            E[n] = (flow_.mixtures().thermalTable().ea_rho(rho[n], T[n], flow_.mixtures().Yi(), n) +
                    real(0.5) * v_[n].magSqr());
        }
    }
    HurMPI::allReduce(sumRelCorrec, MPI_SUM);
    pseudoTimes().cfl().setRelativeCorrection(
        sqrt(sumRelCorrec / real(mesh().allMeshCellNumber())));

    HurMPI::allReduce(countNan, MPI_SUM);
    if (countNan != 0) {
        if (countNan > 50) {
            errorAbortStr(("Too many NaN: " + toString(countNan)));
        } else {
            checkWarningStr(("NaN occurs in results: " + toString(countNan)));
        }
    }
}

void OpenHurricane::LUSGS::fluxJacoDQ(const integer i, const vector &n, realArray &adq) {
    if (i >= mesh().nCells() && dq_[i][0] == real(0) && dq_[i][1] == real(0)) {
        adq = real(0.0);
        return;
    }
    const realArray &dqi = dq_[i];
    solver_.Acdq(i, n, dqi, adq);
}

void OpenHurricane::LUSGS::marching() {
    lowerLoop();

    upperLoop();
    update();
}

void OpenHurricane::LUSGS::lowerLoop(cellRealArray &_cellQ, const vector2DArray &_ra,
                                     const realArray &_dt, realArray &_dq) {
    const auto &_mesh = _cellQ.mesh();
    const auto &_fcl = _mesh.faces();
    const auto &_cell = _mesh.cells();
    const auto nCells = _mesh.nCells();
    const auto &_meshVol = _mesh.cellVolume();
    for (integer n = 0; n < nCells; ++n) {
        real dm = _meshVol[n] / _dt[n] + _cellQ.diagSource()[n];
        real fluxq = Zero;
        for (integer fi = 0; fi < _cell[n].faceSize(); fi++) {
            const integer faceI = _cell[n].facei(fi);
            const auto &cl = _fcl[faceI].leftCell();
            const auto &cr = _fcl[faceI].rightCell();

            const integer m = cl + cr - n;

            if (m < n) {
                const real sl = real(abs(cl - n));
                const real sr = real(abs(cr - n));
                const real coek = _ra[faceI][0] * sr / max(real(1.0), sr) +
                                  _ra[faceI][1] * sl / max(real(1.0), sl);
                fluxq += coek * _dq[m];
            }
        }

        _dq[n] = (_cellQ.rhs()[n] - fluxq) / dm;
    }
}

void OpenHurricane::LUSGS::upperLoop(cellRealArray &_cellQ, const vector2DArray &_ra,
                                     const realArray &_dt, realArray &_dq) {
    const auto &_mesh = _cellQ.mesh();
    const auto &_fcl = _mesh.faces();
    const auto &_cell = _mesh.cells();
    const auto nCells = _mesh.nCells();
    const auto &_meshVol = _mesh.cellVolume();
    for (integer n = nCells - 1; n >= 0; --n) {
        real dm = _meshVol[n] / _dt[n] + _cellQ.diagSource()[n];
        real fluxq = Zero;
        for (integer fi = 0; fi < _cell[n].faceSize(); fi++) {
            const integer faceI = _cell[n].facei(fi);
            const auto &cl = _fcl[faceI].leftCell();
            const auto &cr = _fcl[faceI].rightCell();

            const integer m = cl + cr - n;

            if (m > n) {
                const real sl = real(abs(cl - n));
                const real sr = real(abs(cr - n));
                const real coek = _ra[faceI][0] * sr / max(real(1.0), sr) +
                                  _ra[faceI][1] * sl / max(real(1.0), sl);
                fluxq += coek * _dq[m];
            }
        }

        _dq[n] -= fluxq / dm;
    }
}
