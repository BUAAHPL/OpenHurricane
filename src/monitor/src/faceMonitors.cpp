/*!
 * \file faceMonitors.cpp
 * \brief Main subroutines for monitoring face.
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

#include "faceMonitors.hpp"
#include "controllerSwitch.hpp"
#include "fVArraysInclude.hpp"

namespace OpenHurricane {
    createClassName(faceMonitors);
    registerObjFty(monitor, faceMonitors, controller);
} // namespace OpenHurricane

void OpenHurricane::faceMonitors::readFromCont(const controller &resCont) {
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

void OpenHurricane::faceMonitors::getZoneId() {
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

OpenHurricane::real OpenHurricane::faceMonitors::getFaceMassFlowRate() const {
    const auto &faces = mesh().faces();
    const auto &fW = mesh().faceWeight();
    const auto &fA = mesh().faceArea();
    const auto &fzl = mesh().faceZones();
    const auto &rho = mesh().findObject<cellRealArray>("rho");
    const auto &v = mesh().findObject<cellVectorArray>("v");
    real massFlowRate = Zero;
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const auto fzi = zoneIdList_[i];

        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            const auto wl = fW[fi];

            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();

            const auto rhol = rho[cl] * wl + rho[cr] * (1.0 - wl);
            const auto vl = v[cl] * wl + v[cr] * (real(1.0) - wl);

            massFlowRate += rhol * (vl * fA[fi]);
        }
    }
    return massFlowRate;
}

void OpenHurricane::faceMonitors::writeAndPrintReal(const cellRealArray &phi) const {
    auto &fos = const_cast<fileOsstream &>(fosResidual_);
    if (perFaceZone_) {
        const auto &fzl = mesh().faceZones();
        std::cout.setf(std::ios::showpoint);
        for (integer i = 0; i < zoneIdList_.size(); ++i) {
            const auto si = surIntPtr_[i].surIntegral(phi);
            if (HurMPI::master()) {
                if (printToScreen_) {
                    string name = monitorName_ + "(" + fzl[zoneIdList_[i]].name() + ")";
                    std::cout << "    " << std::left << std::setfill(' ') << std::setw(12)
                              << name.c_str() << ": ";
                    std::cout << std::left << std::setfill(' ') << std::setw(12)
                              << std::setprecision(5) << si;
                    std::cout << "    (" << surIntPtr_[i].printInfo();
                    if (dependOnVar_) {
                        std::cout << " of " << phi.name().c_str();
                    }
                    std::cout << ")" << std::endl;
                }
                if (writeToFile_) {
                    fos.os() << '\t' << std::setprecision(5) << si;
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
        const auto si = surIntPtr_[0].surIntegral(phi);
        if (HurMPI::master()) {
            if (printToScreen_) {
                std::cout << "    " << std::left << std::setfill(' ') << std::setw(12)
                          << monitorName_.c_str() << ": ";
                std::cout << std::left << std::setfill(' ') << std::setw(12) << std::setprecision(5)
                          << si;
                std::cout << "    (" << surIntPtr_[0].printInfo();
                if (dependOnVar_) {
                    std::cout << " of " << phi.name().c_str();
                }
                std::cout << ")" << std::endl;
            }
            if (writeToFile_) {
                fos.os() << '\t' << std::setprecision(5) << si;
                fos.os() << std::endl;
            }
        }
        std::cout.unsetf(std::ios::showpoint);
        std::cout << std::right;
        return;
    }
}

void OpenHurricane::faceMonitors::writeAndPrintVector(const cellVectorArray &phi) const {
    auto &fos = const_cast<fileOsstream &>(fosResidual_);
    if (perFaceZone_) {
        const auto &fzl = mesh().faceZones();
        std::cout.setf(std::ios::showpoint);
        for (integer i = 0; i < zoneIdList_.size(); ++i) {
            const auto si = surIntPtr_[i].surIntegral(phi);
            if (HurMPI::master()) {
                if (printToScreen_) {
                    string name = monitorName_ + "(" + fzl[zoneIdList_[i]].name() + ")";
                    std::cout << "    " << std::left << std::setfill(' ') << std::setw(12)
                              << name.c_str() << ": (";
                    std::cout << std::left << std::setfill(' ') << std::setw(12)
                              << std::setprecision(5) << si[0] << ", ";
                    std::cout << std::left << std::setfill(' ') << std::setw(12)
                              << std::setprecision(5) << si[1] << ", ";
                    std::cout << std::left << std::setfill(' ') << std::setw(12)
                              << std::setprecision(5) << si[2];
                    std::cout << ")    (" << surIntPtr_[i].printInfo();
                    if (dependOnVar_) {
                        std::cout << " of " << phi.name().c_str();
                    }
                    std::cout << ")" << std::endl;
                }
                if (writeToFile_) {
                    fos.os() << '\t' << std::setprecision(5) << si[0] << '\t'
                             << std::setprecision(5) << si[1] << '\t' << std::setprecision(5)
                             << si[2];
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
        const auto si = surIntPtr_[0].surIntegral(phi);
        if (HurMPI::master()) {
            if (printToScreen_) {
                std::cout << "    " << std::left << std::setfill(' ') << std::setw(12)
                          << monitorName_.c_str() << ": (";
                std::cout << std::left << std::setfill(' ') << std::setw(12) << std::setprecision(5)
                          << si[0] << ", ";
                std::cout << std::left << std::setfill(' ') << std::setw(12) << std::setprecision(5)
                          << si[1] << ", ";
                std::cout << std::left << std::setfill(' ') << std::setw(12) << std::setprecision(5)
                          << si[2];
                std::cout << ")    (" << surIntPtr_[0].printInfo();
                if (dependOnVar_) {
                    std::cout << " of " << phi.name().c_str();
                }
                std::cout << ")" << std::endl;
            }
            if (writeToFile_) {
                fos.os() << '\t' << std::setprecision(5) << si[0] << '\t' << std::setprecision(5)
                         << si[1] << '\t' << std::setprecision(5) << si[2];
                fos.os() << std::endl;
            }
        }
        std::cout.unsetf(std::ios::showpoint);
        std::cout << std::right;
        return;
    }
}

void OpenHurricane::faceMonitors::writeAndPrint(const integer step) const {
    auto &fos = const_cast<fileOsstream &>(fosResidual_);
    if (HurMPI::master()) {
        if (writeToFile_) {
            fos.os() << step;
        }
    }

    if (!dependOnVar_) {
        writeAndPrintReal(cellRealArray::nullObject());
    } else {
        if (!mesh().foundOnlyObject(varName_)) {
            errorAbortStr(("Cannot find " + varName_));
        }

        const auto &obj = mesh().findOnlyObject(varName_);
        if (obj.nElements() == 1) {
            const cellRealArray &phi = mesh().findObject<cellRealArray>(varName_);
            writeAndPrintReal(phi);
        } else if (obj.nElements() == 3) {
            const cellVectorArray &phi = mesh().findObject<cellVectorArray>(varName_);
            writeAndPrintVector(phi);
        } else {
            errorAbortStr(
                ("Type with " + toString(obj.nElements()) + " components not supported yet"));
        }
    }
}

OpenHurricane::faceMonitors::faceMonitors(const iteration &iter, const runtimeMesh &mesh,
                                          const controller &cont, const string &name)
    : monitor(iter, mesh, cont, name), zoneNameList_(), zoneIdList_(), surIntPtr_(),
      perFaceZone_(false), dependOnVar_(false),
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
        /*string pw = optCont.findWord("perFaceZone");
        trim(pw);
        stringToUpperCase(pw);
        if (pw == "ON")
        {
                perFaceZone_ = true;
        }
        else if (pw == "OFF")
        {
                perFaceZone_ = false;
        }
        else
        {
                errorAbortStr("Unknown type: " + pw + " in " + optCont.name(), HUR_FUNCTION);
        }*/
    }

    if (perFaceZone_) {
        surIntPtr_.resize(zoneIdList_.size());
        for (integer i = 0; i < zoneIdList_.size(); ++i) {
            integerList zd(integer(1), zoneIdList_[i]);
            surIntPtr_.set(i, surfaceIntegrals::creator(iter, mesh, cont, zd).release());
        }
    } else {
        surIntPtr_.append(surfaceIntegrals::creator(iter, mesh, cont, zoneIdList_).release());
    }
    if (zoneIdList_.size() != 0) {
        dependOnVar_ = surIntPtr_[0].dependOnVariables();
    }

    if (dependOnVar_) {
        if (cont.found("var")) {
            varName_ = cont.findWord("var");
        } else {
            LFatal("var must be given in %s", cont.name().c_str());
        }
    }
}

void OpenHurricane::faceMonitors::monitoring() const {
    if (iter().cStep() % updateStep_ == 0) {
        writeAndPrint(iter().totalStep());
    }
}

void OpenHurricane::faceMonitors::subMonitoring() const {
    if (iter().subIter().cSubStep() % updateStep_ == 0) {
        writeAndPrint(iter().subIter().totalStep());
    }
}