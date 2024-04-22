/*!
 * \file pointMonitors.cpp
 * \brief Main subroutines for monitoring points.
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

#include "pointMonitors.hpp"
#include "cellArrays.hpp"

namespace OpenHurricane {
    createClassName(pointMonitors);
    registerObjFty(monitor, pointMonitors, controller);
} // namespace OpenHurricane

void OpenHurricane::pointMonitors::sendToMaster(realArray &data) const {
    if (HurMPI::master() && HurMPI::isThisProc(processId_)) {
        return;
    }
    if (componentSize_ == 0 || data.size() == 0) {
        return;
    }
    if (HurMPI::master()) {
        if (data.size() != componentSize_) {
            data.resize(componentSize_);
        }
    }

    if (HurMPI::isThisProc(processId_)) {
        if (data.size() != componentSize_) {
            errorAbortStr(("The components size: " + toString(componentSize_) +
                           " is not equal to the size of data array: " + toString(data.size())));
        }
        HurMPI::sendList(data, HurMPI::masterNo(), 1, HurMPI::getComm());
    } else if (HurMPI::master()) {
        HurMPI::Status status;
        HurMPI::recvList(data, processId_, 1, HurMPI::getComm(), &status);
    }
}

void OpenHurricane::pointMonitors::writeOutScreenInMaster(const realArray &data) const {
    if (HurMPI::master()) {
        std::cout << std::endl;
        std::cout << std::left << std::setfill(' ') << std::setw(15);
        std::cout << std::left << std::setfill(' ') << std::setw(10) << monitorName_.c_str();
        std::cout.setf(std::ios::showpoint);
        integer end = componentId_[0];
        integer start;
        for (integer i = 0; i < monitorVarName_.size(); i++) {
            std::cout << std::left << std::setfill(' ') << std::setw(12)
                      << monitorVarName_[i].c_str() << ": ";
            start = end;
            end = componentId_[i + 1];
            if ((end - start) == 1) {
                std::cout << std::left << std::setfill(' ') << std::setw(12) << std::setprecision(5)
                          << data[start];
            } else {
                std::cout << "( ";
                for (integer j = start; j < end; ++j) {
                    std::cout << std::left << std::setfill(' ') << std::setw(12)
                              << std::setprecision(5) << data[j];
                    if (j != end - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << ") ";
            }
        }
        std::cout.unsetf(std::ios::showpoint);
        std::cout << std::right;
    }
}

void OpenHurricane::pointMonitors::writeOutFileInMaster(const realArray &data) const {
    if (HurMPI::master()) {
        fosPoint_.os() << iter().totalStep();
        if (iter().hasPhysicalTimeStep()) {
            fosPoint_.os() << '\t' << toString(iter().pTStep().totalTime()).c_str();
        }
        for (integer i = 0; i < data.size(); ++i) {
            fosPoint_.os() << '\t' << toString(data[i]).c_str();
        }
        fosPoint_.os() << std::endl;
    }
}

void OpenHurricane::pointMonitors::setMonitorVarCmpt() {
    if (monitorVarName_.size() == 0) {
        return;
    }

    integer m = 0;
    componentId_[0] = 0;
    for (integer i = 0; i < monitorVarName_.size(); ++i) {
        if (mesh().foundOnlyObject(monitorVarName_[i])) {
            const auto &obj = mesh().findOnlyObject(monitorVarName_[i]);
            m += obj.nElements();
            componentId_[i + 1] = m;
            monitorVarCmptName_ += obj.outputVarName();
        }
        //else if(...) For vector or tensor variable
        // ...
        //
        else {
            m++;
            componentId_[i + 1] = m;
            monitorVarCmptName_ += monitorVarName_[i];
        }
    }
    componentSize_ = m;
}

void OpenHurricane::pointMonitors::getMonitorVarCmpt(realArray &data) const {
    // Only in the process containing the monitored point.
    if (HurMPI::isThisProc(processId_)) {
        if (monitorVarName_.size() == 0) {
            return;
        }
        integer end = componentId_[0];
        integer start;
        for (integer i = 0; i < monitorVarName_.size(); ++i) {
            start = end;
            end = componentId_[i + 1];

            if (mesh().foundOnlyObject(monitorVarName_[i])) {
                if ((end - start) == 1) {
                    const auto &var = mesh().findObject<cellRealArray>(monitorVarName_[i]);
                    data[start] = var[cellIndex_];
                } else if ((end - start) == 3) {
                    const auto &var = mesh().findObject<cellVectorArray>(monitorVarName_[i]);
                    data[start] = var[cellIndex_][0];
                    data[start + 1] = var[cellIndex_][1];
                    data[start + 2] = var[cellIndex_][2];
                }
            }
            //else if(...) For vector or tensor variable
            // ...
            //
            else if (monitorVarName_[i] == "Ma") {
                data[start] = pointMa();
            }
        }
    }
}

OpenHurricane::real OpenHurricane::pointMonitors::pointMa() const {
    return real();
}

OpenHurricane::pointMonitors::pointMonitors(const iteration &iter, const runtimeMesh &mesh,
                                            const controller &cont, const string &name)
    : monitor(iter, mesh, cont, name), cellIndex_(-1),
      position_(cont.findOrDefault<vector>("position", vector(0))), processId_(), monitorVarName_(),
      monitorVarCmptName_(), componentSize_(0), componentId_(),
      fosPoint_(IOsstream::streamFormat::ASCII_FORMAT, IOsstream::openOptions::ONLY_MASTER,
                std::ios_base::app)

{
    cellIndex_ = mesh.findPointCell(position_);

    integerList idl(HurMPI::getProcSize(), -1);
    processId_ = HurMPI::getProcRank();

    idl[HurMPI::getProcRank()] = cellIndex_;
    HurMPI::allGatherList(idl);
    bool isSet = false;
    for (integer ip = 0; ip < idl.size(); ++ip) {
        if (idl[ip] >= 0) {
            if (!isSet) {
                if (!HurMPI::isThisProc(ip)) {
                    cellIndex_ = idl[ip];
                    processId_ = ip;
                }
                isSet = true;
            }
        }
    }
    if (cellIndex_ == -1) {
        errorAbortStr(("The position: " + toString(position_) + " is not within the cell domain"));
    }

    std::string rnl = cont.findText("varList");
    replaceAllMarks(rnl, "\n", " ");
    if (!rnl.empty()) {
        size_t pos = 0;
        stdStringList rll;
        split(rnl, rll, ",");
        monitorVarName_.resize(rll.size());
        for (integer i = 0; i < rll.size(); ++i) {
            monitorVarName_[i] = trimCopy(rll[i]);
        }
    }
    removeDuplicate(monitorVarName_);
    componentId_.resize(monitorVarName_.size() + 1, -1);
    setMonitorVarCmpt();
    if (writeToFile_) {
        fosPoint_.changeFileName(outFile_);
        if (!iter.restart() || (iter.restart() && !fosPoint_.existed())) {
            fosPoint_.changeMode(std::ios_base::out);
        }
        fosPoint_.open();
        if (!iter.restart() || (iter.restart() && !fosPoint_.existed())) {
            if (iter.isSteadyFlow()) {
                fosPoint_.os() << "variables = iter," << monitorVarCmptName_.c_str() << std::endl;
            } else {
                fosPoint_.os() << "variables = iter,time_s," << monitorVarCmptName_.c_str()
                               << std::endl;
            }
        }
    } else {
        fosPoint_.close();
    }
}

void OpenHurricane::pointMonitors::monitoring() const {
    realArray data(componentSize_);
    getMonitorVarCmpt(data);
    sendToMaster(data);

    if (printToScreen_) {
        writeOutScreenInMaster(data);
    }
    if (writeToFile_) {
        writeOutFileInMaster(data);
    }
}