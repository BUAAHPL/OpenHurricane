/*!
 * \file residuals.cpp
 * \brief Main subroutines for monitoring residuals.
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

#include "residuals.hpp"
#include "controllerSwitch.hpp"
#include "fVArraysInclude.hpp"

namespace OpenHurricane {
    createClassName(residuals);
    registerObjFty(monitor, residuals, controller);
} // namespace OpenHurricane

void OpenHurricane::residuals::readFromCont(const controller &resCont) {
    if (updateStep_ < 1) {
        checkWarningStr(("Update step:" + toString(updateStep_) +
                         " cannot less than 1, and has been set to 1"));
    }
    if (!resCont.found("varList")) {
        setDefaultMonitorVars();
        return;
    }
    residualsNameList_.clear();
    std::string rnl = resCont.findText("varList");
    replaceAllMarks(rnl, "\n", " ");
    if (!rnl.empty()) {
        size_t pos = 0;
        stdStringList rll;
        split(rnl, rll, ",");
        residualsNameList_.resize(rll.size());
        for (integer i = 0; i < rll.size(); ++i) {
            residualsNameList_[i] = trimCopy(rll[i]);
        }
    }
    removeDuplicate(residualsNameList_);
}

void OpenHurricane::residuals::setDefaultMonitorVars() {
    residualsNameList_.append("rho");
    residualsNameList_.append("u");
    residualsNameList_.append("v");
    residualsNameList_.append("w");
    residualsNameList_.append("E");
}

void OpenHurricane::residuals::parsingList(const controller &resCont, const integerList &cmptNum,
                                           stringList &cmptName) {
    integer totalNum = 0;
    for (integer cmpti = 0; cmpti < cmptNum.size(); cmpti++) {
        totalNum += cmptNum[cmpti];
    }

    if (resCont.found("residualConvergence")) {
        checkResidualConvergence_ =
            controllerSwitch(resCont)("residualConvergence", checkResidualConvergence_);
        if (checkResidualConvergence_) {
            convergenceThreshold_.resize(totalNum, real(0));
            integer m = 0;
            for (integer i = 0; i < residualsNameList_.size(); ++i) {
                if (cmptNum[i] == 1) {
                    convergenceThreshold_[m] =
                        resCont.findOrDefault<real>(residualsNameList_[i], 1e-6);
                    m++;
                } else if (cmptNum[i] == 3) {
                    vector var = resCont.findOrDefault<vector>(residualsNameList_[i], vector(1e-6));
                    convergenceThreshold_[m++] = var.x();
                    convergenceThreshold_[m++] = var.y();
                    convergenceThreshold_[m++] = var.z();
                } else {
                    errorAbortStr(("Invalid component number: " + toString(cmptNum[i])));
                }
            }
        } else {
            convergenceThreshold_.clear();
        }
    }

    integer m = 0;
    stringList newName(totalNum);
    for (integer i = 0; i < residualsNameList_.size(); ++i) {
        if (cmptNum[i] == 1) {
            replaceAllMarks(cmptName[i], "\"", "");
            newName[m] = cmptName[i];
            m++;
        } else if (cmptNum[i] == 3) {
            std::string str = cmptName[i];
            replaceAllMarks(str, "\"", " ");
            if (!str.empty()) {
                size_t pos = 0;
                stdStringList rll;
                split(str, rll, ",");
                for (integer j = 0; j < rll.size(); ++j) {
                    newName[m++] = trimCopy(rll[j]);
                }
            }
        } else {
            errorAbortStr(("Invalid component number: " + toString(cmptNum[i])));
        }
    }
    cmptName.transfer(newName);
}

OpenHurricane::residuals::residuals(const iteration &iter, const runtimeMesh &mesh,
                                    const controller &cont, const string &name)
    : monitor(iter, mesh, cont, name), residualsNameList_(), checkResidualConvergence_(false),
      convergenceThreshold_(),
      fosResidual_(IOsstream::streamFormat::ASCII_FORMAT, IOsstream::openOptions::ONLY_MASTER,
                   std::ios_base::out) {
    readFromCont(cont);
    resCmptNum_ = mesh.initResidualsListAndCheck(residualsNameList_, residualsCmptNameList_);
    parsingList(cont, resCmptNum_, residualsCmptNameList_);

    if (writeToFile_) {
        fosResidual_.changeFileName(outFile_);

        fosResidual_.open();
    } else {
        fosResidual_.close();
    }
}

bool OpenHurricane::residuals::resConvergence() const {
    if (!checkResidualConvergence_) {
        return false;
    }

    bool flag = true;
    if (HurMPI::master()) {
        integer m = 0;
        for (integer i = 0; i < resCmptNum_.size(); i++) {
            if (resCmptNum_[i] == 1) {
                const real curRhs =
                    mesh().findObject<cellRealArray>(residualsNameList_[i]).curRhsAve();
                if (curRhs > convergenceThreshold_[m]) {
                    flag = false;
                    break;
                }
                m++;
            } else if (resCmptNum_[i] == 3) {
                const auto curRhs =
                    mesh().findObject<cellVectorArray>(residualsNameList_[i]).curRhsAve();
                if (curRhs.x() > convergenceThreshold_[m++]) {
                    flag = false;
                    break;
                }
                if (curRhs.y() > convergenceThreshold_[m++]) {
                    flag = false;
                    break;
                }
                if (curRhs.z() > convergenceThreshold_[m++]) {
                    flag = false;
                    break;
                }
            } else {
                errorAbortStr(("Invalid component number: " + toString(resCmptNum_[i])));
            }
        }
    }
    HurMPI::bcast(&flag, 1, feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
    return flag;
}

void OpenHurricane::residuals::writeResiduals(fileOsstream &fos) const {
    const geometryMesh &mesh = iter().findObject<geometryMesh>(iter().name());

    if (iter().cStep() == updateStep_ || ((iter().cStep() - 1) % 100 == 0)) {
        if (HurMPI::master()) {
            std::cout << std::endl;
            std::cout << std::left << std::setfill(' ') << std::setw(10) << "Iter";
            for (integer i = 0; i < residualsCmptNameList_.size(); i++) {
                std::cout << std::left << std::setfill(' ') << std::setw(12)
                          << residualsCmptNameList_[i].c_str();
            }
            std::cout << std::left << std::setfill(' ') << std::setw(12) << "Remaining";
            if (iter().hasPhysicalTimeStep() && !iter().hasSubIteration()) {
                std::cout << std::left << std::setfill(' ') << std::setw(20) << "Current Time(s)";
                std::cout << std::left << std::setfill(' ') << std::setw(20) << "Total Time(s)";
                std::cout << std::left << std::setfill(' ') << std::setw(20) << "Time step(s)";
            }
            std::cout << std::endl;
        }
    }
    if (HurMPI::master()) {
        // Write the residual results to screen.
        std::cout << std::left << std::setfill(' ') << std::setw(10) << iter().totalStep();
        if (fos.opened()) {
            fos.os() << iter().totalStep();
        }
    }

    mesh.writeResiduals(fos, updateStep_);

    if (HurMPI::master()) {
        std::cout << std::left << std::setfill(' ') << std::setw(12)
                  << iter().maxStep() - iter().cStep();
        if (iter().hasPhysicalTimeStep() && !iter().hasSubIteration()) {
            std::cout.setf(std::ios::showpoint);
            std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(5)
                      << iter().pTStep().totalTime();
            std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(5)
                      << iter().pTStep().maxTime();
            std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(5)
                      << iter().pTStep().pTimeStep();
            std::cout.unsetf(std::ios::showpoint);
        }
        std::cout << std::endl;
        std::cout << std::right;
        if (fos.opened()) {
            fos.os() << std::right;
            fos.os() << std::endl;
            fos.flush();
        }
    }
    if (!iter().hasSubIteration()) {
        if (resConvergence()) {
            iter().setResConveg();
        }
    }
}

void OpenHurricane::residuals::writeSubIterResiduals(fileOsstream &fos) const {
    if (!iter().hasSubIteration()) {
        return;
    }

    const geometryMesh &mesh = iter().findObject<geometryMesh>(iter().name());

    if (iter().subIter().cSubStep() == updateStep_) {
        if (HurMPI::master()) {
            std::cout << "   " << std::left << std::setfill(' ') << std::setw(16) << "Sub-iter";
            for (integer i = 0; i < residualsCmptNameList_.size(); i++) {
                std::cout << std::left << std::setfill(' ') << std::setw(12)
                          << residualsCmptNameList_[i].c_str();
            }
            std::cout << std::left << std::setfill(' ') << std::setw(12) << "Remaining"
                      << std::endl;
        }
    }
    if (HurMPI::master()) {
        std::cout << "   " << std::left << std::setfill(' ') << std::setw(16)
                  << iter().subIter().cSubStep();
        if (fos.opened()) {
            fos.os() << "   " << iter().subIter().totalStep();
        }
    }
    mesh.writeResiduals(fos, updateStep_);

    if (HurMPI::master()) {
        std::cout << std::left << std::setfill(' ') << std::setw(16)
                  << iter().subIter().maxSubStep() - iter().subIter().cSubStep() << std::endl;
        std::cout << std::right;
        if (fos.opened()) {
            fos.os() << std::endl;
            fos.flush();
        }
    }

    if (resConvergence()) {
        iter().subIter().setResConveg();
    }
}
