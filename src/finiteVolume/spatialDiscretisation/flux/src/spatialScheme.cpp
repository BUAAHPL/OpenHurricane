/*!
 * \file spatialScheme.cpp
 * \brief Main subroutines for spatial scheme.
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

#include "spatialScheme.hpp"

namespace OpenHurricane {
    createClassNameStr(spatialScheme, "spatialScheme");
    createObjFty(spatialScheme, controller);
} // namespace OpenHurricane

OpenHurricane::spatialScheme::spatialScheme(const controller &cont, const runtimeMesh &mesh,
                                            flowModel &flow)
    : mesh_(mesh), flows_(flow), rhoFlux_(object("rhoFlux", mesh, object::NOT_WRITE), mesh),
      v_(flow.v()), thermo_(flow.thermo()), reconstrPtr_(nullptr), objectList_(), countParam_(0),
      paramMap_(mesh.primitiveParamSize()), flux_(5, Zero), shockFactor_(flow.shockFactor()),
      isLowMachCorr_(false) {
    isLowMachCorr_ = controllerSwitch(cont)("isLowMachCorrection", isLowMachCorr_);

    reconstrPtr_ = reconstruction::creator(cont);
}

OpenHurricane::uniquePtr<OpenHurricane::spatialScheme>
OpenHurricane::spatialScheme::creator(const controller &cont, const runtimeMesh &mesh,
                                      flowModel &flow) {
    string spschemeType = cont.findWord(spatialScheme::className_);

    Pout << "    Info: setting spatial scheme: " << spschemeType << std::endl;
    defineInObjCreator(spatialScheme, static_cast<std::string>(spschemeType), controller,
                       (cont, mesh, flow));
}

OpenHurricane::uniquePtr<OpenHurricane::spatialScheme>
OpenHurricane::spatialScheme::creator(const string &schemeType, const controller &cont,
                                      const runtimeMesh &mesh, flowModel &flow) {
    Pout << "    Info: setting spatial scheme: " << schemeType << std::endl;
    defineInObjCreator(spatialScheme, static_cast<std::string>(schemeType), controller,
                       (cont, mesh, flow));
}

OpenHurricane::spatialScheme::~spatialScheme() noexcept {
    for (integer i = 0; i < objectList_.size(); ++i) {
        objectList_[i] = nullptr;
    }
    reconstrPtr_.clear();
}

void OpenHurricane::spatialScheme::basicFlux() {
    const auto &mesh = v_.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();

    const auto &fA = mesh.faceArea();
    const auto &fC = mesh.faceCntr();
    const auto &cC = mesh.cellCentre();

    auto &v = v_;
    auto &rho = const_cast<rhoThermo &>(thermo_).rho();
    auto &p = const_cast<rhoThermo &>(thermo_).p();
    auto &T = const_cast<rhoThermo &>(thermo_).T();
    auto &E = const_cast<rhoThermo &>(thermo_).E();
    auto &gama = const_cast<rhoThermo &>(thermo_).gamma();

    rho.rhs() = Zero;
    v.rhs() = Zero;
    T.rhs() = Zero;
    E.rhs() = Zero;

    reconstrPtr_->calcGrad(rho);
    reconstrPtr_->calcGrad(v);
    reconstrPtr_->calcGrad(p);
    reconstrPtr_->calcGrad(T);
    reconstrPtr_->calcGrad(E);
    reconstrPtr_->calcGrad(gama);

    reconstrPtr_->calcLimiter(rho);
    reconstrPtr_->calcLimiter(v);
    reconstrPtr_->calcLimiter(p);
    reconstrPtr_->calcLimiter(E);
    reconstrPtr_->calcLimiter(gama);

    mixture &mixtures = const_cast<rhoThermo &>(thermo_).mixtures();

    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        // Symmetric plane
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const integer &cl = fl[fi].leftCell();

                const vector vFlux = p[cl] * fA[fi];
                v.rhs()[cl] += vFlux;
                rhoFlux_[fi] = Zero;
            }
        }
        //Other cell face inviscid flux calculations
        else {
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                real rhol, rhor;
                reconstrPtr_->calcReconstruction(rho, fi, rhol, rhor);
                real pl, pr;
                vector vl, vr;
                real gaml, gamr;
                real el, er;
                reconstrPtr_->calcReconstruction(p, fi, pl, pr);
                if (rhol <= 0.0 || rhor <= 0.0 || pl <= 0.0 || pr <= 0.0) {
                    pl = p[cl];
                    pr = p[cr];
                    rhol = rho[cl];
                    rhor = rho[cr];
                    vl = v[cl];
                    vr = v[cr];
                    gaml = gama[cl];
                    gamr = gama[cr];
                    el = E[cl];
                    er = E[cr];
                } else {
                    reconstrPtr_->calcReconstruction(v, fi, vl, vr);
                    reconstrPtr_->calcReconstruction(gama, fi, gaml, gamr);
                    reconstrPtr_->calcReconstruction(E, fi, el, er);
                    if (isLowMachCorr_) {
                        lowMachCorrection(rhol, rhor, gaml, gamr, vl, vr, pl, pr, fA[fi]);
                    }
                }

                const real al = sqrt(gaml * pl / rhol);
                const real ar = sqrt(gamr * pr / rhor);
                flux_ = Zero;
                real blend = min(shockFactor_[cl], shockFactor_[cr]);

                calcFlux(rhor, rhol, vr, vl, pr, pl, gamr, gaml,
                         //er + pr / rhor - 0.5 * vr.magSqr(), el + pl / rhol - 0.5 * vl.magSqr(),
                         er, el, ar, al, fA[fi], blend, flux_);

                rhoFlux_[fi] = flux_[0];
                if (std::isnan(flux_[0]) || isinf(flux_[0])) {
                    std::string errMsg;
                    errMsg = "density inv_flux is not a number\n";
                    errMsg += " rhol = " + toString(rhol);
                    errMsg += " rhor = " + toString(rhor) + "\n";
                    errMsg += " pl = " + toString(pl);
                    errMsg += " pr = " + toString(pr) + "\n";
                    errMsg += " In face center: " + toString(fC[fi]);
                    errorAbortStr(errMsg);
                }
                if (std::isnan(flux_[4]) || isinf(flux_[4])) {
                    std::string errMsg;
                    errMsg = "energy inv_flux is not a number\n";
                    errMsg += " el = " + toString(el);
                    errMsg += " er = " + toString(er) + "\n";
                    errMsg += " In face center: " + toString(fC[fi]);
                    errorAbortStr(errMsg);
                }
                rho.rhs()[cl] += flux_[0];
                v.rhs()[cl][0] += flux_[1];
                v.rhs()[cl][1] += flux_[2];
                v.rhs()[cl][2] += flux_[3];
                E.rhs()[cl] += flux_[4];

                if (cr < mesh.nCells()) {
                    rho.rhs()[cr] -= flux_[0];
                    v.rhs()[cr][0] -= flux_[1];
                    v.rhs()[cr][1] -= flux_[2];
                    v.rhs()[cr][2] -= flux_[3];
                    E.rhs()[cr] -= flux_[4];
                }
            }
        }
    }
}

void OpenHurricane::spatialScheme::invFlux(cellRealArray &cellQ) const {
    invFluxTemplate(cellQ);
}

void OpenHurricane::spatialScheme::invFlux(cellRealArray &cellQ, const realBounded &bound) const {
    cellQ.rhs() = Zero;
    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);
    const auto &mesh = cellQ.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();
    const auto &fC = mesh.faceCntr();

    const auto &fA = mesh.faceArea();
    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            continue;
        }
        //Other cell face inviscid flux calculations
        else {
            realArray ql(fzl[fzi].size());
            realArray qr(fzl[fzi].size());
            reconstrPtr_->calcReconstruction(cellQ, fzi, ql, qr);
            integer count = 0;
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();

                real flux;
                if (rhoFlux_[fi] >= 0.0) {
                    flux = rhoFlux_[fi] * bound.bounding(qr[count], cellQ[cr]);
                } else {
                    flux = rhoFlux_[fi] * bound.bounding(ql[count], cellQ[cl]);
                }
                //Add the numerical flux to the right-hand-term of both cells sharing the face
                cellQ.rhs()[cl] += flux;
                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= flux;
                }
                count++;
            }
        }
    }
}

void OpenHurricane::spatialScheme::invFlux(cellVectorArray &cellQ,
                                           const vectorBounded &bound) const {
    cellQ.rhs() = Zero;
    reconstrPtr_->calcGrad(cellQ);
    reconstrPtr_->calcLimiter(cellQ);
    const auto &mesh = cellQ.mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();
    const auto &fC = mesh.faceCntr();

    const auto &fA = mesh.faceArea();
    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            continue;
        }
        //Other cell face inviscid flux calculations
        else {
            vectorArray ql(fzl[fzi].size());
            vectorArray qr(fzl[fzi].size());
            reconstrPtr_->calcReconstruction(cellQ, fzi, ql, qr);
            integer count = 0;
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();

                vector flux;
                if (rhoFlux_[fi] >= 0.0) {
                    flux = rhoFlux_[fi] * bound.bounding(qr[count], cellQ[cr]);
                } else {
                    flux = rhoFlux_[fi] * bound.bounding(ql[count], cellQ[cl]);
                }
                //Add the numerical flux to the right-hand-term of both cells sharing the face
                cellQ.rhs()[cl] += flux;
                if (cr < mesh.nCells()) {
                    cellQ.rhs()[cr] -= flux;
                }
                count++;
            }
        }
    }
}

void OpenHurricane::spatialScheme::invFlux(cellVectorArray &cellQ) const {
    invFluxTemplate(cellQ);
}

void OpenHurricane::spatialScheme::invFluxSpecies(PtrList<cellRealArray> &yi,
                                                  const bool withLastSpc) const {
    if (withLastSpc) {
        for (integer i = 0; i < yi.size(); ++i) {
            yi[i].rhs() = Zero;
            reconstrPtr_->calcGrad(yi[i]);
            reconstrPtr_->calcLimiter(yi[i]);
        }
    } else {
        for (integer i = 0; i < yi.size() - 1; ++i) {
            yi[i].rhs() = Zero;
            reconstrPtr_->calcGrad(yi[i]);
            reconstrPtr_->calcLimiter(yi[i]);
        }
    }
    const auto &mesh = yi[0].mesh();
    const auto &fzl = mesh.faceZones();
    const auto &fl = mesh.faces();
    const auto &fC = mesh.faceCntr();

    const auto &fA = mesh.faceArea();
    for (integer fzi = 0; fzi < fzl.size(); fzi++) {
        //Wall cell faces invicid fluxes: only explicit pressure fluxes considered.
        if (fzl[fzi].isWall() || fzl[fzi].isSymmetric()) {
            continue;
        }
        //Other cell face inviscid flux calculations
        else {
            realArrayArray ql;
            realArrayArray qr;
            if (withLastSpc) {
                ql.resize(yi.size());
                qr.resize(yi.size());
                for (integer i = 0; i < yi.size(); ++i) {
                    ql[i].resize(fzl[fzi].size());
                    qr[i].resize(fzl[fzi].size());
                    reconstrPtr_->calcReconstruction(yi[i], fzi, ql[i], qr[i]);
                }
            } else {
                ql.resize(yi.size() - 1);
                qr.resize(yi.size() - 1);
                for (integer i = 0; i < yi.size() - 1; ++i) {
                    ql[i].resize(fzl[fzi].size());
                    qr[i].resize(fzl[fzi].size());
                    reconstrPtr_->calcReconstruction(yi[i], fzi, ql[i], qr[i]);
                }
            }

            integer count = 0;
            for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); fi++) {
                const auto &cl = fl[fi].leftCell();
                const auto &cr = fl[fi].rightCell();
                if (withLastSpc) {
                    real yilt = Zero;
                    real yirt = Zero;
                    for (integer i = 0; i < yi.size(); ++i) {
                        ql[i][count] = min(real(1), max(real(0), ql[i][count]));
                        qr[i][count] = min(real(1), max(real(0), qr[i][count]));
                        yilt += ql[i][count];
                        yirt += qr[i][count];
                    }
                    for (integer i = 0; i < yi.size(); ++i) {
                        ql[i][count] /= max(tiny, yilt);
                        qr[i][count] /= max(tiny, yirt);
                        real flux;
                        if (rhoFlux_[fi] >= 0.0) {
                            flux = rhoFlux_[fi] * qr[i][count];
                        } else {
                            flux = rhoFlux_[fi] * ql[i][count];
                        }
                        //Add the numerical flux to the right-hand-term of both cells sharing the face
                        yi[i].rhs()[cl] += flux;
                        if (cr < mesh.nCells()) {
                            yi[i].rhs()[cr] -= flux;
                        }
                    }
                } else {
                    real yilt = Zero;
                    real yirt = Zero;
                    for (integer i = 0; i < yi.size() - 1; ++i) {
                        ql[i][count] = min(real(1), max(real(0), ql[i][count]));
                        qr[i][count] = min(real(1), max(real(0), qr[i][count]));
                        yilt += ql[i][count];
                        yirt += qr[i][count];
                    }
                    real yiLastL = 1.0 - yilt;
                    real yiLastR = 1.0 - yirt;

                    for (integer i = 0; i < yi.size() - 1; ++i) {
                        if (yiLastL < 0) {
                            ql[i][count] /= yilt;
                        }
                        if (yiLastR < 0) {
                            qr[i][count] /= yirt;
                        }
                        real flux;
                        if (rhoFlux_[fi] >= 0.0) {
                            flux = rhoFlux_[fi] * qr[i][count];
                        } else {
                            flux = rhoFlux_[fi] * ql[i][count];
                        }
                        //Add the numerical flux to the right-hand-term of both cells sharing the face
                        yi[i].rhs()[cl] += flux;
                        if (cr < mesh.nCells()) {
                            yi[i].rhs()[cr] -= flux;
                        }
                    }
                }
                count++;
            }
        }
    }
}
