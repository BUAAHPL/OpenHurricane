/*!
 * \file freeStream.cpp
 * \brief The subroutines and functions of basic free streams
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

#include "freeStream.hpp"
#include "flowModel.hpp"

void OpenHurricane::freeStream::autoSetFreeStream(const controller &bbCont) {
#ifdef HUR_DEBUG
    Pout << "    Info: setting free stream automatically..." << std::endl;
#endif // HUR_DEBUG
    bool isSet = false;
    auto &mapEle = bbCont.mapEntries();
    for (auto &e : mapEle) {
        const auto &bcCont = bbCont.subController(e.first);
        const auto bcKey = bcCont.findWord("bcType");
        //=====================================================
        //   Pressure far field
        //=====================================================
        if (bcKey == "pressureFarField" || bcKey == "pressure-far-field" ||
            bcKey == "supersonicInlet") {
#ifdef HUR_DEBUG
            Pout << "    Info: setting free stream with boundary: " << e.first << std::endl;
#endif // HUR_DEBUG
            setFreeStream(bcCont);
            isSet = true;
            break;
        } else if (bcKey == "pressureInlet" || bcKey == "pressure-inlet" ||
                   bcKey == "subsonicInlet") {
#ifdef HUR_DEBUG
            Pout << "    Info: setting free stream with boundary: " << e.first << std::endl;
#endif // HUR_DEBUG
            setFreeStream(bcCont);
            isSet = true;
            break;
        } else if (bcKey == "massFlowInlet" || bcKey == "mass-flow-inlet") {
#ifdef HUR_DEBUG
            Pout << "    Info: setting free stream with boundary: " << e.first << std::endl;
#endif // HUR_DEBUG
            setFreeStream(bcCont);
            isSet = true;
            break;
        }
    }
    if (!isSet) {
        LFatal("The free stream cannot be set automatically, please check!");
    }
}

void OpenHurricane::freeStream::setFreeStream(const controller &bcCont) {
    getSpeciesMassFractions(bcCont);

    if (!bcCont.found("p")) {
        LFatal("Cannot find pressure for setting free stream in %s", bcCont.name().c_str());
    }
    if (!bcCont.found("T")) {
        LFatal("Cannot find temperature for setting free stream in %s", bcCont.name().c_str());
    }
    p_ = bcCont.findType<real>("p", p_);
    T_ = bcCont.findType<real>("T", T_);

    rho_ = flow_.mixtures().thermalTable().eos().rhom(p_, T_, yi_);

    gamma_ = flow_.mixtures().thermalTable().gamma(p_, T_, yi_);
    mu_ = flow_.mixtures().transTable().mu(p_, T_, yi_);

    if (bcCont.found("v")) {
        vMag_ = bcCont.findType<real>("v", vMag_);

        Ma_ = vMag_ / sqrt(gamma_ * p_ / rho_);
    } else if (bcCont.found("ma")) {
        Ma_ = bcCont.findType<real>("ma", Ma_);
        vMag_ = Ma_ * sqrt(gamma_ * p_ / rho_);
    } else {
        std::string errMsg;
        errMsg = "No velocity specification defined in freeStream: ";
        errMsg += bcCont.name();
        errMsg += ".\nPlease check freeStream.";
        errorAbortStr(errMsg);
    }
    if (bcCont.found("direction")) {
        const std::regex contFlag("(.*?)\\s*\\((.*?)(?:,\\s*?)(.*?)(?:,\\s*?)(.*?)\\s*?\\)");
        std::smatch what;
        string str = bcCont.findWord("direction");
        std::regex_search(str, what, contFlag);
        std::string type = what[1];
        std::istringstream cmpt1(what[2]);
        std::istringstream cmpt2(what[3]);
        std::istringstream cmpt3(what[4]);
        real x = static_cast<real>(feature<real>(cmpt1));
        real y = static_cast<real>(feature<real>(cmpt2));
        real z = static_cast<real>(feature<real>(cmpt3));
        vector direction(x, y, z);

        if (type == "car") /*!Cartesian coordinates.*/
        {
            v_ = direction.normalized() * vMag_;
        }
    } else {
        vector direction(1.0, 0, 0);
        v_ = direction * vMag_;
    }
}

void OpenHurricane::freeStream::getSpeciesMassFractions(const controller &bcCont) {
    const auto &species = flow_.mixtures().species();
    yi_.resize(species.size(), Zero);

    if (species.size() == 1) {
        yi_ = 1.0;
    } else {
        const controller &specCont = bcCont.subController("species");
        integer givenFlag = 0; // 0 = mass fraction; 1 = mole fraction.
        if (specCont.found("givenBy")) {
            const auto &flagW = specCont.findWord("givenBy");
            if (flagW == "massFraction") {
                givenFlag = 0;
            } else if (flagW == "moleFraction") {
                givenFlag = 1;
            } else {
                errorAbortStr(
                    ("Unknown type: " + flagW + " for given species in: " + specCont.name()));
            }
        }

        real yisum = Zero;
        for (integer i = 0; i < species.size() - 1; ++i) {
            std::string spn = species[i].name();
            yi_[i] = specCont.findOrDefault<real>(spn, 0.0);
            if (yi_[i] < 0) {
                errorAbortStr(("Value must not lower than zero for species: " + spn + " in " +
                               specCont.name()));
            }
            yisum += yi_[i];
        }
        if (yisum > 1.0) {
            LFatal("The summation of all species must not greater than one in %s",
                   specCont.name().c_str());
        }

        yi_[species.size() - 1] = 1.0 - yisum;

        if (givenFlag == 1) {
            species.Xi2Yi(yi_, yi_);
        }
    }
}

OpenHurricane::freeStream::freeStream(const controller &cont, const flowModel &flow)
    : flow_(flow), Ma_(0.0), v_(Zero), vMag_(0.0), p_(0.0), T_(0.0), rho_(0.0), mu_(0.0),
      gamma_(0.0), yi_(1) {
    bool autoSet = false;
    const auto &topCont = cont.topCont();
    if (topCont.found("freeStream")) {
        const auto &fsCont = topCont.subController("freeStream");
        if (fsCont.found("specifiedType")) {
            string sptw = fsCont.findWord("specifiedType");
            if (sptw == "auto") {
                autoSet = true;
            } else if (sptw == "givenByFaceZone") {
                string bcName = fsCont.subController("givenByFaceZone").findWord("givenBy");
                const auto &bbCont = topCont.subController("boundaryCondition");
                setFreeStream(bbCont.subController(bcName));
            } else if (sptw == "givenByParameters") {
                const auto &bbCont = fsCont.subController("givenByParameters");
                setFreeStream(bbCont);
            } else {
                errorAbortStr(("Unknown options: " + sptw + " in " + fsCont.name()));
            }
        } else {
            autoSet = true;
        }
    } else {
        autoSet = true;
    }

    if (autoSet) {
        if (topCont.found("boundaryCondition")) {
            const auto &bbCont = topCont.subController("boundaryCondition");
            autoSetFreeStream(bbCont);
        } else {
            LFatal("Cannot set free stream because the boundary condition cannot be found in %s",
                   topCont.name().c_str());
        }
    }
}