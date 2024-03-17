/*!
 * \file writeAveFieldVar.cpp
 * \brief Main subroutines for writing average field variables.
 * \author Chen Zhenyi
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

#include "writeAveFieldVar.hpp"
#include "calculateFieldVar.hpp"

namespace OpenHurricane {
    createClassNameStr(writeAveFieldVar, "writeAveFieldVar");
    registerObjFty(writeFieldVar, writeAveFieldVar, controller);
} // namespace OpenHurricane

void OpenHurricane::writeAveFieldVar::changeOutVarNameToAve() {
    for (integer i = 0; i < writeFieldVarList_.size(); ++i) {
        writeFieldVarList_[i] = writeFieldVarList_[i] + "_Ave";
    }
}

OpenHurricane::writeAveFieldVar::writeAveFieldVar(const flowModel &flows, const iteration &iter,
                                                  const controller &cont, const string &writeId,
                                                  std::map<std::string, object *> &outFieldVarMap)
    : writeFieldVar(flows, iter, cont, writeId, outFieldVarMap), timeSumField(), startTime_(0),
      curretTimeGap_(0), endTimeAveCal_(0) {
    changeOutVarNameToAve();
    timeSumField::setTimeSumVarList(flows, mesh(), writeFieldVarList_);
    setoutVarMap(writeFieldVarList_);
    if (iter_.hasPhysicalTimeStep()) {
        startTime_ = iter_.pTStep().totalTime();
    }
}

void OpenHurricane::writeAveFieldVar::updating() {
    if (!writeNow()) {
        return;
    }
    if (!iter_.hasPhysicalTimeStep()) {
        return;
    }
    setCurretTimeGap();
    for (integer i = 0; i < writeFieldVarList_.size(); ++i) {
        const auto &nam = writeFieldVarList_[i];
        auto pfw = getPrefix(nam);
        bool isSpecies = false;
        for (integer j = 0; j < yi().size(); ++j) {
            if (pfw == yi()[j].name()) {
                calcTimeAveSpecies();
                setOutputField(yi()[j].name() + "_Ave", yiAve()[j]);
                isSpecies = true;
                break;
            }
        }
        if (isSpecies) {
            continue;
        }
        if (nam == "Ma_Ave") {
            setOutputField(
                nam, calculateFieldVar::calcMachNumber(calcTimeAveVar(v()), calcTimeAveVar(gama()),
                                                       calcTimeAveVar(p()), calcTimeAveVar(rho())));
        } else if (nam == "pt_Ave") {
            calcTotalPressure("Ma_Ave", "pt_Ave", calcTimeAveVar(p()), calcTimeAveVar(gama()),
                              calcTimeAveVar(v()), calcTimeAveVar(rho()));

        } else if (nam == "Tt_Ave") {
            if (flows_.mixtures().isSingular()) {
                calcTotalTemperature("Ma_Ave", "Tt_Ave", calcTimeAveVar(T()),
                                     calcTimeAveVar(gama()), calcTimeAveVar(v()),
                                     calcTimeAveVar(p()), calcTimeAveVar(rho()));
            } else {
                calcTimeAveSpecies();

                setOutputField("Tt_Ave", calculateFieldVar::totalTemperature(
                                             flows_, calcTimeAveVar(p()), calcTimeAveVar(T()),
                                             calcTimeAveVar(v()), yiAve()));
            }
        } else {
            for (auto iter2 = mesh().table().begin(); iter2 != mesh().table().end(); ++iter2) {
                if (pfw == iter2->first) {
                    if (iter2->second->nElements() == 1) {
                        cellRealArray *f = static_cast<cellRealArray *>(iter2->second);
                        setOutputField(nam, calcTimeAveVar(*f));
                    }
                    if (iter2->second->nElements() == 3) {
                        cellVectorArray *f = static_cast<cellVectorArray *>(iter2->second);
                        setOutputField(nam, calcTimeAveVar(*f));
                    }
                    break;
                }
            }
        }
    }
}

void OpenHurricane::writeAveFieldVar::correct() const {
    sumField(mesh());
}

void OpenHurricane::writeAveFieldVar::reseting() const {
    writeFieldVar::reseting();
    reset();
}

void OpenHurricane::writeAveFieldVar::calcTimeAveSpecies() {
    if (yiSet()) {
        return;
    }
    for (integer i = 0; i < yi().size(); ++i) {
        auto oldyi = yiAve().set(
            i, new cellRealArray(object(yi()[i].name() + "_Ave", mesh(), object::TEMPORARY), mesh(),
                                 calcTimeAveVar(yi()[i])));
        HurDelete(oldyi);
    }
}

void OpenHurricane::writeAveFieldVar::setCurretTimeGap() {
    curretTimeGap_ = iter_.pTStep().totalTime() - startTime_;
}

OpenHurricane::fileName OpenHurricane::writeAveFieldVar::getFileName() const {
    fileName outN = iter_.outputName();
    auto pathOut = outN.parentPath();
    auto fname = outN.name(true);
    string fext;
    fext = ".plt";

    fname += "-";
    fname += "Average";
    fname += "-";
    fname += toString(iter_.totalStep());
    outN = fname + fext;
    outN = pathOut / outN;

    return outN;
}

bool OpenHurricane::writeAveFieldVar::writeNow() const {
    if (!isWritting_ || writeFieldVarList_.size() == 0 || iter_.isSteadyFlow()) {
        return false;
    } else {
        if (iter_.cStep() % iter_.writeOutputStep() == 0 ||
            (iter_.end() && !argParse::noWriteAtEnd())) {
            return true;
        }
    }
    return false;
}
