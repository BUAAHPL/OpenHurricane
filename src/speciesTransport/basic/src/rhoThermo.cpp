/*!
 * \file rhoThermo.cpp
 * \brief Main subroutines for the rho thermo.
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
#include "rhoThermo.hpp"
#include "getBoundariesFromController.hpp"

OpenHurricane::rhoThermo::rhoThermo(const runtimeMesh &mesh)
    : p_(object("p", mesh, object::WRITE_RELAY_OUTPUT), mesh),
      T_(object("T", mesh, object::WRITE_RELAY_OUTPUT), mesh),
      rho_(object("rho", mesh, object::WRITE_RELAY_OUTPUT, object::PRIMITIVE), mesh),
      E_(object("E", mesh, object::WRITE_RELAY_OUTPUT, object::PRIMITIVE), mesh),
      gamma_(object("gamma", mesh, object::WRITE_OUTPUT), mesh),
      cp_(object("cp", mesh, object::WRITE_OUTPUT), mesh), mixPtr_(nullptr) {
    mixPtr_.reset(new mixture(mesh, 1));
}

OpenHurricane::rhoThermo::rhoThermo(const runtimeMesh &mesh, const controller &cont,
                                const bool inviscous)
    : p_(object("p", mesh, object::WRITE_RELAY_OUTPUT), mesh),
      T_(object("T", mesh, object::WRITE_RELAY_OUTPUT), mesh),
      rho_(object("rho", mesh, object::WRITE_RELAY_OUTPUT, object::PRIMITIVE), mesh),
      E_(object("E", mesh, object::WRITE_RELAY_OUTPUT, object::PRIMITIVE), mesh),
      gamma_(object("gamma", mesh, object::WRITE_OUTPUT), mesh),
      cp_(object("cp", mesh, object::WRITE_OUTPUT), mesh), mixPtr_(nullptr) {
    mixPtr_.reset(new mixture(mesh, cont.subController("mixture")));
    getBoundariesFromController::setBoundariesController(const_cast<controller &>(cont.topCont()),
                                                         *mixPtr_);
}

void OpenHurricane::rhoThermo::THa(const cellRealArray &ha, const cellRealArray &rho,
                               cellRealArray &T) const {
    const integer nCell = T.mesh().nCells();

    for (integer cellI = 0; cellI < nCell; ++cellI) {
        integer iflag = -1;
        T[cellI] = mixPtr_->thermalTable().THa_rho(ha[cellI], rho[cellI], T[cellI], iflag,
                                                   mixPtr_->Yic(cellI));

        if (iflag == 0 || iflag == 2) {
            const_cast<cellRealArray &>(ha)[cellI] =
                mixPtr_->thermalTable().ha_rho(rho_[cellI], T[cellI], mixPtr_->Yic(cellI));
        }
    }
}

void OpenHurricane::rhoThermo::THa(const cellRealArray &ha, const cellRealArray &rho, cellRealArray &T,
                               integerArray &_flag) const {
    const integer nCell = T.mesh().nCells();

    for (integer cellI = 0; cellI < nCell; ++cellI) {
        _flag[cellI] = -1;
        T[cellI] = mixPtr_->thermalTable().THa_rho(ha[cellI], rho[cellI], T[cellI], _flag[cellI],
                                                   mixPtr_->Yic(cellI));

        if (_flag[cellI] == 0 || _flag[cellI] == 2) {
            const_cast<cellRealArray &>(ha)[cellI] =
                mixPtr_->thermalTable().ha_rho(rho_[cellI], T[cellI], mixPtr_->Yic(cellI));
        }
    }
}

void OpenHurricane::rhoThermo::TEa(const cellRealArray &ea, const cellRealArray &rho,
                               cellRealArray &T) const {
    const integer nCell = T.mesh().nCells();

    for (integer cellI = 0; cellI < nCell; ++cellI) {
        integer iflag = -1;
        T[cellI] = mixPtr_->thermalTable().TEa_rho(ea[cellI], rho[cellI], T[cellI], iflag,
                                                   mixPtr_->Yic(cellI));

        if (iflag == 0 || iflag == 2) {
            const_cast<cellRealArray &>(ea)[cellI] =
                mixPtr_->thermalTable().ea_rho(rho_[cellI], T[cellI], mixPtr_->Yic(cellI));
        }
    }
}

void OpenHurricane::rhoThermo::TEa(const cellRealArray &ea, const cellRealArray &rho, cellRealArray &T,
                               integerArray &_flag) const {
    const integer nCell = T.mesh().nCells();

    for (integer cellI = 0; cellI < nCell; ++cellI) {
        _flag[cellI] = -1;
        T[cellI] = mixPtr_->thermalTable().TEa_rho(ea[cellI], rho[cellI], T[cellI], _flag[cellI],
                                                   mixPtr_->Yic(cellI));

        if (_flag[cellI] == 0 || _flag[cellI] == 2) {
            const_cast<cellRealArray &>(ea)[cellI] =
                mixPtr_->thermalTable().ea_rho(rho_[cellI], T[cellI], mixPtr_->Yic(cellI));
        }
    }
}

void OpenHurricane::rhoThermo::THa(const cellRealArray &Ha, const cellVectorArray &v,
                               const cellRealArray &rho, cellRealArray &T) const {
    const integer nCell = T.mesh().nCells();

    for (integer cellI = 0; cellI < nCell; ++cellI) {
        integer iflag = -1;
        T[cellI] =
            mixPtr_->thermalTable().THa_rho((Ha[cellI] - real(0.5) * v[cellI].magSqr()), rho[cellI],
                                            T[cellI], iflag, mixPtr_->Yic(cellI));

        if (iflag == 0 || iflag == 2) {
            const_cast<cellRealArray &>(Ha)[cellI] =
                (mixPtr_->thermalTable().ha_rho(rho_[cellI], T[cellI], mixPtr_->Yic(cellI)) +
                 real(0.5) * v[cellI].magSqr());
        }
    }
}

void OpenHurricane::rhoThermo::THa(const cellRealArray &Ha, const cellVectorArray &v,
                               const cellRealArray &rho, cellRealArray &T,
                               integerArray &_flag) const {
    const integer nCell = T.mesh().nCells();

    for (integer cellI = 0; cellI < nCell; ++cellI) {
        _flag[cellI] = 1;
        T[cellI] =
            mixPtr_->thermalTable().THa_rho((Ha[cellI] - real(0.5) * v[cellI].magSqr()), rho[cellI],
                                            T[cellI], _flag[cellI], mixPtr_->Yic(cellI));

        if (_flag[cellI] == 0 || _flag[cellI] == 2) {
            const_cast<cellRealArray &>(Ha)[cellI] =
                (mixPtr_->thermalTable().ha_rho(rho_[cellI], T[cellI], mixPtr_->Yic(cellI)) +
                 real(0.5) * v[cellI].magSqr());
        }
    }
}

void OpenHurricane::rhoThermo::TEa(const cellRealArray &Ea, const cellVectorArray &v,
                               const cellRealArray &rho, cellRealArray &T) const {
    const integer nCell = T.mesh().nCells();

    for (integer cellI = 0; cellI < nCell; ++cellI) {
        integer iflag = -1;

        T[cellI] =
            mixPtr_->thermalTable().TEa_rho((Ea[cellI] - real(0.5) * v[cellI].magSqr()), rho[cellI],
                                            T[cellI], iflag, mixPtr_->Yi(), cellI);

        if (iflag == 0 || iflag == 2) {
            const_cast<cellRealArray &>(Ea)[cellI] =
                (mixPtr_->thermalTable().ea_rho(rho_[cellI], T[cellI], mixPtr_->Yi(), cellI) +
                 real(0.5) * v[cellI].magSqr());
        }
    }
}

void OpenHurricane::rhoThermo::TEa(const cellRealArray &Ea, const cellVectorArray &v,
                               const cellRealArray &rho, cellRealArray &T,
                               integerArray &_flag) const {
    const integer nCell = T.mesh().nCells();

    for (integer cellI = 0; cellI < nCell; ++cellI) {
        _flag[cellI] = 1;
        real rhoe = Ea[cellI] * rho[cellI];
        real ek = 0.5 * v[cellI].magSqr();
        T[cellI] = mixPtr_->thermalTable().TEa_rho((Ea[cellI] - ek), rho[cellI], T[cellI],
                                                   _flag[cellI], mixPtr_->Yi(), cellI);
        real Td = T[cellI];
        if (_flag[cellI] == 0 || _flag[cellI] == 2) {

            const_cast<cellRealArray &>(Ea)[cellI] =
                (mixPtr_->thermalTable().ea_rho(rho_[cellI], T[cellI], mixPtr_->Yi(), cellI) + ek);
        }
    }
}

void OpenHurricane::rhoThermo::correctLimits(cellRealArray &T, const integer cellI, const real lowLmt,
                                         const real highLmt, const integer flag) const {
    const integer nCell = T.mesh().nCells();
    const auto &cells = T.mesh().cells();
    const auto &faces = T.mesh().faces();
    const auto &cellCtr = T.mesh().cellCentre();
    if (cellI >= nCell || cellI < 0) {
        return;
    }

    if (flag == 1) {
        return;
    } else if (flag == 0) {
        real Tsum = 0;
        integer count = 0;
        real twgt = 0;
        for (integer fi = 0; fi < cells[cellI].faceSize(); fi++) {
            const integer fj = cells[cellI].facei(fi);
            const integer cl = faces[fj].leftCell();
            const integer cr = faces[fj].rightCell();
            const integer m = cl + cr - cellI;
            if (m < nCell) {
                if (T[m] > lowLmt) {
                    const real wgt = 1.0 / dist(cellCtr[m], cellCtr[cellI]);
                    Tsum += wgt * T[m];
                    twgt += wgt;
                    count++;
                }
            }
        }
        if (count > 0) {
            T[cellI] = Tsum / twgt;
        } else {
            T[cellI] = lowLmt;
        }
    } else if (flag == 2) {
        real Tsum = 0;
        integer count = 0;
        real twgt = 0;
        for (integer fi = 0; fi < cells[cellI].faceSize(); fi++) {
            const integer fj = cells[cellI].facei(fi);
            const integer cl = faces[fj].leftCell();
            const integer cr = faces[fj].rightCell();
            const integer m = cl + cr - cellI;
            if (m < nCell) {
                if (T[m] < highLmt) {
                    const real wgt = 1.0 / dist(cellCtr[m], cellCtr[cellI]);
                    Tsum += wgt * T[m];
                    twgt += wgt;
                    count++;
                }
            }
        }
        if (count > 0) {
            T[cellI] = Tsum / twgt;
        } else {
            T[cellI] = highLmt;
        }
    }
}

void OpenHurricane::rhoThermo::p(const cellRealArray &rho, const cellRealArray &T,
                             cellRealArray &p) const {
    integer nCell = T.mesh().nCells();

    for (integer cellI = 0; cellI < nCell; ++cellI) {
        p[cellI] = mixPtr_->thermalTable().eos().pm(rho[cellI], T[cellI], mixPtr_->Yi(), cellI);
    }
}

void OpenHurricane::rhoThermo::mu(const cellRealArray &p, const cellRealArray &T,
                              cellRealArray &mu) const {
    integer nCell = T.mesh().nCells();

    for (integer cellI = 0; cellI < nCell; ++cellI) {
        mu[cellI] = mixPtr_->transTable().mu(p[cellI], T[cellI], mixPtr_->Yic(cellI));
    }
}