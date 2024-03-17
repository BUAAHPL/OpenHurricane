/*!
 * \file RANSModel.cpp
 * \brief Main subroutines for the RANS turbulence model.
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

#include "RANSModel.hpp"

void OpenHurricane::RANSModel::limitTurbVarIncrease(cellRealArray &cellQ, realArray &oldQ) const {
    for (integer n = 0; n < mesh().nCells(); ++n) {
        cellQ[n] = max(min(cellQ[n], (real(1.0) + coefTurb_) * oldQ[n]),
                       (real(1.0) - coefTurb_) * oldQ[n]);
    }
}

OpenHurricane::RANSModel::RANSModel(const controller &cont, flowModel &ev)
    : turbulenceModel(cont, ev), intensity_(10), viscosityRatio_(0), lengthScale_(0),
      hydraulicD_(0), coefTurb_(cont.findOrDefault<real>("coefTurb", 0.05)),
      bndSpMethod_(intensityAndLength), averageLimit_(false) {
    turbulenceModel::yPtr_.reset(new wallDistance(ev.mesh(), cont));
}

void OpenHurricane::RANSModel::inletBndSetting(controller &turbCont) {
    if (turbCont.findWord("specificationMethod") == "intensityAndLength") {
        bndSpMethod_ = intensityAndLength;
        if (turbCont.found("intensity")) {
            intensity_ = turbCont.findType<real>("intensity", intensity_) * real(0.01);
            if (turbCont.found("length-scale")) {
                lengthScale_ = turbCont.findType<real>("length-scale", lengthScale_);
            } else {
                LFatal(
                    "Turbulent parameters specification found in boundary type: %s are not enough.",
                    turbCont.findWord("bcType").c_str());
            }
        } else {
            LFatal("No turbulent parameter specification found in boundary type: %s",
                   turbCont.findWord("bcType").c_str());
        }
    } else if (turbCont.findWord("specificationMethod") == "intensityAndHydraulicD") {
        bndSpMethod_ = intensityAndHYdraulicD;
        if (turbCont.found("intensity")) {
            intensity_ = turbCont.findType<real>("intensity", intensity_) * real(0.01);
            if (turbCont.found("Hydraulic-Diameter")) {
                hydraulicD_ = turbCont.findType<real>("Hydraulic-Diameter", hydraulicD_);
            } else {
                LFatal(
                    "Turbulent parameters specification found in boundary type: %s are not enough.",
                    turbCont.findWord("bcType").c_str());
            }
        } else {
            LFatal("No turbulent parameter specification found in boundary type: %s",
                   turbCont.findWord("bcType").c_str());
        }
    } else if (turbCont.findWord("specificationMethod") == "viscosityRatio") {
        bndSpMethod_ = viscosityRatio;
        if (turbCont.found("viscosity-ratio")) {
            viscosityRatio_ = turbCont.findType<real>("viscosity-ratio", viscosityRatio_);
            if (turbCont.found("viscosity-ratio") &&
                cont().subController("flow")
                        .subController("turbulence")
                        .findWord("turbulenceModel") == "SpalartAllmaras") {
                LFatal("Turbulent parameters specification found in boundary type: %s do not "
                       "support this method in SpalartAllmaras model.",
                       turbCont.findWord("bcType").c_str());
            } else if (turbCont.found("intensity")) {
                intensity_ = turbCont.findType<real>("intensity", intensity_) * real(0.01);
            } else {
                LFatal(
                    "Turbulent parameters specification found in boundary type: %s are not enough.",
                    turbCont.findWord("bcType").c_str());
            }
        } else {
            LFatal("No turbulent parameter specification found in boundary type: %s",
                   turbCont.findWord("bcType").c_str());
        }
    } else if (turbCont.findWord("specificationMethod") == "origianalTurbEquation") {
        bndSpMethod_ = origianalTurbEquation;
    } else if (turbCont.findWord("specificationMethod") == "givenDirectly") {
        bndSpMethod_ = givenDirectly;
    } else {
        LFatal("Turbulent parameters specification found in boundary type: %s are not enough.",
               turbCont.findWord("bcType").c_str());
    }
}

void OpenHurricane::RANSModel::limitNegative(cellRealArray &cellQ, const real lowestV,
                                         const bool average, const bool reportToScreen) {
    const auto &fL = mesh().faces();
    const auto &cs = mesh().cells();
    const auto &cc = mesh().cellCntr();
    const cellIntegerArray &CFLFlag = mut().tb().findObject<cellIntegerArray>("CFLFlag");
    if (lowestV < 0.0) {
        LFatal("Negative limit is not allow: %e", lowestV);
    }
    integer countExc = 0;
    for (integer n = 0; n < mesh().nCells(); ++n) {
        if (cellQ[n] < lowestV) {
            if (CFLFlag[n] == 1) {
                const_cast<cellIntegerArray &>(CFLFlag)[n] = 0;
            }
            countExc++;
            cellQ[n] = lowestV;
            if (average) {
                integer count = 0;
                real sumQ = Zero;
                real sumW = Zero;
                for (integer i = 0; i < cs[n].faceSize(); ++i) {
                    const auto fi = cs[n].facei(i);

                    const auto cl = fL[fi].leftCell();
                    const auto cr = fL[fi].rightCell();

                    const integer m = cl + cr - n;

                    if (m < mesh().nCells()) {
                        if (cellQ[m] > lowestV) {
                            const real ft = 1.0 / dist(cc[m], cc[n]);
                            sumQ += ft * cellQ[m];
                            sumW += ft;
                            count++;
                        }
                    }
                }
                if (count > 0) {
                    cellQ[n] = sumQ / sumW;
                }
            }
        }
    }
    HurMPI::allReduce(countExc, MPI_SUM);
    if (countExc > 0 && reportToScreen) {
        if (HurMPI::master()) {
            std::cerr << "    Variable: " << cellQ.name() << " limited to " << std::scientific
                      << lowestV << " in " << countExc << " cells" << std::endl;
        }
    }
}
