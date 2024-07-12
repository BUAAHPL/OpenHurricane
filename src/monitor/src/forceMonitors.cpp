/*!
 * \file forceMonitors.cpp
 * \brief Main subroutines for monitoring force.
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

#include "forceMonitors.hpp"
#include "controllerSwitch.hpp"
#include "fVArraysInclude.hpp"

namespace OpenHurricane {
    createClassName(forceMonitors);
    registerObjFty(monitor, forceMonitors, controller);
} // namespace OpenHurricane

void OpenHurricane::forceMonitors::readFromCont(const controller &resCont) {
    if (updateStep_ < 1) {
        checkWarningStr(("Update step:" + toString(updateStep_) +
                         " cannot less than 1, and has been set to 1"));
    }
    if (!resCont.found("zoneList")) {
        LFatal("Cannot find monitored faces in %s", resCont.name().c_str());
    }
    zoneNameList_.clear();
    std::string rnl = resCont.findText("zoneList");
    replaceAllMarks(rnl, "\n", " ");
    if (!rnl.empty()) {
        size_t pos = 0;
        stdStringList rll;
        split(rnl, rll, ",");
        zoneNameList_.resize(rll.size());
        for (integer i = 0; i < rll.size(); ++i) {
            zoneNameList_[i] = trimCopy(rll[i]);
        }
    }
    removeDuplicate(zoneNameList_);
    getZoneId();
}

void OpenHurricane::forceMonitors::getZoneId() {
    const auto &fzl = mesh().faceZones();
    zoneIdList_.resize(zoneNameList_.size());
    integer count = 0;
    for (integer i = 0; i < zoneNameList_.size(); ++i) {
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].name() == zoneNameList_[i]) {
                zoneIdList_[i] = fzi;
                count++;
                break;
            } else if (fzi == fzl.size() - 1) {
                checkWarningStr((zoneNameList_[i] + " not found in face zone list"));
            }
        }
    }

    zoneIdList_.resize(count);
}

void OpenHurricane::forceMonitors::writeAndPrintReal() const {
    auto &fos = const_cast<fileOsstream &>(fosResidual_);
    if (perFaceZone_) {
        const auto &fzl = mesh().faceZones();
        std::cout.setf(std::ios::showpoint);
        for (integer i = 0; i < zoneIdList_.size(); ++i) {
            forcesPtr_[i].computing();
            const auto si = forcesPtr_[i].reportArray();
            if (HurMPI::master()) {
                if (printToScreen_) {
                    string name = monitorName_ + "(" + fzl[zoneIdList_[i]].name() + ")";
                    std::cout << "    " << std::left << std::setfill(' ') << std::setw(12)
                              << name.c_str() << ": " << std::endl;
                    if (i == 0) {
                        for (integer j = 0; j < forcesPtr_[i].reportName().size(); j++) {
                            if (j == 0) {
                                std::cout << "        ";
                            }
                            std::cout << std::left << std::setfill(' ') << std::setw(12)
                                      << forcesPtr_[i].reportName()[j].c_str();
                        }
                        std::cout << std::endl;
                    }
                    for (integer j = 0; j < si.size(); ++j) {
                        if (j == 0) {
                            std::cout << "        ";
                        }
                        std::cout << std::left << std::setfill(' ') << std::setw(12)
                                  << std::setprecision(5) << si[j];
                    }
                    std::cout << std::endl;
                }
                if (writeToFile_) {
                    for (integer j = 0; j < si.size(); ++j) {
                        fos.os() << '\t' << std::setprecision(5) << si[j];
                    }
                }
            }
        }
        if (HurMPI::master()) {
            if (writeToFile_) {
                fos.os() << std::endl;
            }
        }
        std::cout.unsetf(std::ios::showpoint);
        std::cout << std::right;
        return;
    } else {
        std::cout.setf(std::ios::showpoint);
        forcesPtr_[0].computing();
        const auto si = forcesPtr_[0].reportArray();
        if (HurMPI::master()) {
            if (printToScreen_) {
                string name = monitorName_;
                std::cout << "    " << std::left << std::setfill(' ') << std::setw(12)
                          << name.c_str() << ": " << std::endl;
                for (integer j = 0; j < forcesPtr_[0].reportName().size(); j++) {
                    if (j == 0) {
                        std::cout << "        ";
                    }
                    std::cout << std::left << std::setfill(' ') << std::setw(12)
                              << forcesPtr_[0].reportName()[j].c_str();
                }
                std::cout << std::endl;
                for (integer j = 0; j < si.size(); ++j) {
                    if (j == 0) {
                        std::cout << "        ";
                    }
                    std::cout << std::left << std::setfill(' ') << std::setw(12)
                              << std::setprecision(5) << si[j];
                }
                std::cout << std::endl;
            }
            if (writeToFile_) {
                for (integer j = 0; j < si.size(); ++j) {
                    fos.os() << '\t' << std::setprecision(5) << si[j];
                }
                fos.os() << std::endl;
            }
        }
        std::cout.unsetf(std::ios::showpoint);
        std::cout << std::right;
        return;
    }
}

void OpenHurricane::forceMonitors::writeAndPrint(const integer step) const {
   /* if (iter().hasSubIteration()) {
        if ((iter().subIter().isLooping() && !iter().subIter().isResConvergence()) ||
            iter().subIter().checkMin()) {
            return;
        }
    }*/
    auto &fos = const_cast<fileOsstream &>(fosResidual_);
    if (HurMPI::master()) {
        if (writeToFile_) {
            fos.os() << step;
            if (iter().hasPhysicalTimeStep()) {
                fos.os() << '\t' << toString(iter().pTStep().totalTime()).c_str();
            }
        }
    }
    writeAndPrintReal();
}

OpenHurricane::forceMonitors::forceMonitors(const iteration &iter, const runtimeMesh &mesh,
                                            const controller &cont, const string &name)
    : monitor(iter, mesh, cont, name), zoneNameList_(), zoneIdList_(), perFaceZone_(false),
      fosResidual_(IOsstream::streamFormat::ASCII_FORMAT, IOsstream::openOptions::ONLY_MASTER,
                   std::ios_base::out) {
    readFromCont(cont);

    if (writeToFile_) {
        fosResidual_.changeFileName(outFile_);

        fosResidual_.open();
    } else {
        fosResidual_.close();
    }

    if (cont.found("options")) {
        const auto &optCont = cont.subController("options");
        perFaceZone_ = controllerSwitch(optCont)("perFaceZone", perFaceZone_);
    }

    if (perFaceZone_) {
        forcesPtr_.resize(zoneIdList_.size());
        for (integer i = 0; i < zoneIdList_.size(); ++i) {
            integerList zd(integer(1), zoneIdList_[i]);
            forcesPtr_.set(i, forces::creator(iter, mesh, cont, zd).release());
        }
    } else {
        forcesPtr_.append(forces::creator(iter, mesh, cont, zoneIdList_).release());
    }
}

void OpenHurricane::forceMonitors::monitoring() const {
    if (iter().cStep() % updateStep_ == 0) {
        writeAndPrint(iter().totalStep());
    }
}

void OpenHurricane::forceMonitors::subMonitoring() const {
    // should not monitor in the sub-iteration
    /*if (iter().subIter().cSubStep() % updateStep_ == 0) {
        writeAndPrint(iter().subIter().totalStep());
    }*/
}