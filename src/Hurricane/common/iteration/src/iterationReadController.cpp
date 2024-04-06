/*!
 * \file iterationReadController.cpp
 * \brief The subroutines and functions of CFD time advance iteration
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
#include "geometryMesh.hpp"
#include "hdf5I.hpp"
#include "hdf5O.hpp"
#include "iteration.hpp"
#include "profiles.hpp"
#include "tecplotWriter.hpp"

void OpenHurricane::iteration::readCont() {
    std::string str;
    if (HurMPI::master()) {
        if (report) {
            PLInfo("    Info: reading controller: %s\n", caseFileName().c_str());
        }
        str = readFileToString(caseFileName());
    }
    cont_.readXML(str);
    cont_.add("inputPath", std::string(inputPath()));
    if (!cont_.found("iteration")) {
        LFatal("Cannot find iteration in %s", cont_.name().c_str());
    }

    const auto &interCont = cont_.subController("iteration");
    getIterSteps(interCont);
    getFlowState(interCont);
    getTimeStepType(interCont);
    getPhysicalTimeStep(interCont);
    getSubIterationStep(interCont);
    getRestartFrom(interCont);
    auto &interCont2 = cont_.subController("iteration");
    getProfile(interCont2);
    if (report) {
        PLInfo("    Info: done reading controller: %s\n", caseFileName().c_str());
    }
}

void OpenHurricane::iteration::readCont(std::string contStr) {
    cont_.readXML(contStr);
    cont_.add("inputPath", std::string(inputPath()));
    if (!cont_.found("iteration")) {
        LFatal("Cannot find iteration in %s", cont_.name().c_str());
    }

    const auto &interCont = cont_.subController("iteration");
    getIterSteps(interCont);
    getFlowState(interCont);
    getTimeStepType(interCont);
    getPhysicalTimeStep(interCont);
    getSubIterationStep(interCont);
    getRestartFrom(interCont);
    auto &interCont2 = cont_.subController("iteration");
    getProfile(interCont2);
}

void OpenHurricane::iteration::readAndAddCont(const fileName &contN) {
    auto contFileName = contN;
    if (!contFileName.isAbsolute()) {
        contFileName = inputPath() / contFileName;
    }
    std::string str;
    if (HurMPI::master()) {
        str = readFileToString(caseFileName());
    }
    cont_.readXML(str);
    const auto &interCont = cont_.subController("iteration");
    getIterSteps(interCont);
    getFlowState(interCont);
    getTimeStepType(interCont);
    getPhysicalTimeStep(interCont);
    getSubIterationStep(interCont);
    getRestartFrom(interCont);
    auto &interCont2 = cont_.subController("iteration");
    getProfile(interCont2);
}

void OpenHurricane::iteration::getIterSteps(const controller &interCont) {
    maxStep_ = interCont.findOrDefault<integer>("maxStep", 10000000);
    totalStep_ = interCont.findOrDefault<integer>("totalStep", 0);
    writeOutputStep_ = interCont.findOrDefault<integer>("writeOutputStep", 500);
}

void OpenHurricane::iteration::getFlowState(const controller &interCont) {
    if (interCont.found("flowState")) {
        string flst = interCont.findWord("flowState");
        stringToUpperCase(flst);
        if (flst == "UNSTEADY") {
            flowState_ = flowStateType::unsteady;
        } else if (flst == "STEADY") {
            flowState_ = flowStateType::steady;
        } else {
            LFatal("Unknown flowState: %s in %s", flst.c_str(), interCont.name().c_str());
        }
    }
}

void OpenHurricane::iteration::getTimeStepType(const controller &interCont) {
    if (interCont.found("pseudoTime")) {
        const auto &pseCont = interCont.subController("pseudoTime");
        if (pseCont.found("timeStepType")) {
            string flst = pseCont.findWord("timeStepType");
            if (flst == "globalTimeStep") {
                timestepType_ = timestepType::globalTimeStep;
            } else if (flst == "localTimeStep") {
                timestepType_ = timestepType::localTimeStep;
            } else {
                LFatal("Unknown timeStepType: %s in %s", flst.c_str(), pseCont.name().c_str());
            }
        }
    }
}

void OpenHurricane::iteration::getPhysicalTimeStep(const controller &interCont) {
    if (interCont.found("timeMethod")) {
        const auto &tsCont = interCont.subController("timeMethod");
        if (tsCont.found("physicalTimes")) {
            const auto &phyCont = tsCont.subController("physicalTimes");
            real phyTS = phyCont.findType<real>("phyTimeStep", real(0.0));
            real maxPhyTS = phyCont.findType<real>("maxPhyTime", real(0.0));
            real totalPhyTS = phyCont.findType<real>("totalPhyTime", real(0.0));

            pTStepPtr_.reset(new physicalTimeStep(phyTS, maxPhyTS, totalPhyTS));
            controllerSwitch mycont(phyCont);
            bool whatSet = mycont("dynamicTimeStep", false);
            if (whatSet) {
                pTStepPtr_->setDynamicTimeStep();
                if (phyCont.found("dynamicCFL")) {
                    pTStepPtr_->setDynamicCFL(phyCont.findOrDefault<real>("dynamicCFL", 0.8));
                }
            } else {
                pTStepPtr_->unsetDynamicTimeStep();
            }
        }
    } else if (interCont.found("phyTimeStep")) {
        real phyTS = interCont.findType<real>("phyTimeStep", real(0.0));
        real maxPhyTS = interCont.findType<real>("maxPhyTime", real(0.0));
        real totalPhyTS = interCont.findType<real>("totalPhyTime", real(0.0));

        pTStepPtr_.reset(new physicalTimeStep(phyTS, maxPhyTS, totalPhyTS));
        controllerSwitch mycont(interCont);
        bool whatSet = mycont("dynamicTimeStep", false);
        if (whatSet) {
            pTStepPtr_->setDynamicTimeStep();
            if (interCont.found("dynamicCFL")) {
                pTStepPtr_->setDynamicCFL(interCont.findOrDefault<real>("dynamicCFL", 0.8));
            }
        } else {
            pTStepPtr_->unsetDynamicTimeStep();
        }
    }
}

void OpenHurricane::iteration::getSubIterationStep(const controller &interCont) {
    if (interCont.found("timeMethod")) {
        const auto &tsCont = interCont.subController("timeMethod");
        if (tsCont.found("dualTimeStep")) {
            const auto &dtsCont = tsCont.subController("dualTimeStep");
            if (dtsCont.found("subIteration")) {
                const auto &subIterCont = dtsCont.subController("subIteration");

                integer maxSubTS = subIterCont.findType<integer>("maxSubStep", integer(0));
                integer minSubTS = subIterCont.findType<integer>("minSubStep", integer(0));
                subIterPtr_.reset(new subIteration(maxSubTS, minSubTS));
            }
        }
    } else if (interCont.found("maxSubStep")) {

        integer maxSubTS = interCont.findType<integer>("maxSubStep", integer(0));
        integer minSubTS = interCont.findType<integer>("minSubStep", integer(0));
        subIterPtr_.reset(new subIteration(maxSubTS, minSubTS));
    }
}

void OpenHurricane::iteration::getRestartFrom(const controller &interCont) {
    if (interCont.found("restartFrom")) {
        restart_ = true;
        restartFrom_ = interCont.findWord("restartFrom");
        if (!restartFrom_.isAbsolute()) {
            restartFrom_ = inputPath() / restartFrom_;
        }
    }
}

void OpenHurricane::iteration::getProfile(controller &interCont) {
    if (interCont.found("profile")) {
        auto &prfCont = interCont.subController("profile");
        string prfType = "OpenHurricane";
        if (prfCont.found("profileType")) {
            prfType = prfCont.findWord("profileType");
        }
        if (prfCont.found("profileFile")) {
            fileName profileFile = prfCont.findWord("profileFile");
            if (!profileFile.isAbsolute()) {
                profileFile = inputPath() / profileFile;
            }
            uniquePtr<profiles> myPrfPtr = profiles::creator(profileFile, prfType);
            myPrfPtr->read(prfCont);
        }
    }
}
