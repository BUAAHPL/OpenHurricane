/*!
 * \file writePulseFieldVar.cpp
 * \brief Main subroutines for writing pulse field variables.
 * \author Chen Zhenyi
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

#include "writePulseFieldVar.hpp"
#include "calculateFieldVar.hpp"

namespace OpenHurricane {
    createClassNameStr(writePulseFieldVar, "writePulseFieldVar");
    registerObjFty(writeFieldVar, writePulseFieldVar, controller);
} // namespace OpenHurricane

void OpenHurricane::writePulseFieldVar::changeOutVarNameToPulse() {
    for (integer i = 0; i < writeFieldVarList_.size(); ++i) {
        writeFieldVarList_[i] = writeFieldVarList_[i] + "_Pulse";
    }
}

OpenHurricane::writePulseFieldVar::writePulseFieldVar(
    const flowModel &flows, const iteration &iter, const controller &cont, const string &writeId,
    std::map<std::string, object *> &outFieldVarMap)
    : writeFieldVar(flows, iter, cont, writeId, outFieldVarMap), timeSumField(), startTime_(0),
      curretTimeGap_(0), endTimeAveCal_(0) {
    changeOutVarNameToPulse();
    timeSumField::setTimeSumVarList(flows, mesh(), writeFieldVarList_);
    setoutVarMap(writeFieldVarList_);
    if (iter_.hasPhysicalTimeStep()) {
        startTime_ = iter_.pTStep().totalTime();
    }
}

void OpenHurricane::writePulseFieldVar::updating() {
    if (!writeNow()) {
        return;
    }
    if (!iter_.hasPhysicalTimeStep()) {
        return;
    }
    if (endTimeAveCal_ <= iter_.pTStep().totalTime()) {
        setCurretTimeGap();

        for (integer i = 0; i < writeFieldVarList_.size(); ++i) {
            const auto &nam = writeFieldVarList_[i];
            auto pfw = getPrefix(nam);
            auto Iter1 = mesh().table().find(pfw);
            auto Iter2 = mesh().table().find(pfw + "_Ave");
            if (Iter1 != mesh().table().end() && Iter2 != mesh().table().end()) {
                if (Iter1->second->nElements() == 1) {
                    const cellRealArray *f1 = static_cast<const cellRealArray *>(Iter1->second);
                    const cellRealArray *f2 = static_cast<const cellRealArray *>(Iter2->second);
                    setOutputField(nam, *f2 - (*f1));
                }
                if (Iter1->second->nElements() == 3) {
                    const cellVectorArray *f1 = static_cast<const cellVectorArray *>(Iter1->second);
                    const cellVectorArray *f2 = static_cast<const cellVectorArray *>(Iter2->second);
                    setOutputField(nam, *f2 - (*f1));
                }
            } else if (Iter1 == mesh().table().end() && Iter2 != mesh().table().end()) {
                if (nam == "Ma_Pulse") {
                    const cellRealArray *f2 = static_cast<const cellRealArray *>(Iter2->second);
                    setOutputField(
                        nam, *f2 - calculateFieldVar::calcMachNumber(v(), gama(), p(), rho()));
                }
                if (nam == "pt_Pulse") {
                    const cellRealArray *f2 = static_cast<const cellRealArray *>(Iter2->second);
                    setOutputField(nam, *f2 - calcTotalPressure("Ma", p(), gama(), v(), rho()));
                }
                if (nam == "Tt_Pulse") {
                    const cellRealArray *f2 = static_cast<const cellRealArray *>(Iter2->second);
                    if (flows_.mixtures().isSingular()) {
                        setOutputField(
                            nam, *f2 - calcTotalTemperature("Ma", T(), gama(), v(), p(), rho()));
                    } else {
                        setOutputField(nam, *f2 - calculateFieldVar::totalTemperature(
                                                      flows_, p(), T(), v(), yi()));
                    }
                }
            } else if (Iter1 != mesh().table().end() && Iter2 == p().tb().table().end()) {
                if (nam == "Ma_Pulse") {
                    cellRealArray *f1 = static_cast<cellRealArray *>(Iter1->second);
                    setOutputField(nam, calculateFieldVar::calcMachNumber(
                                            calcTimeAveVar(v()), calcTimeAveVar(gama()),
                                            calcTimeAveVar(p()), calcTimeAveVar(rho())) -
                                            *f1);
                } else if (nam == "pt_Pulse") {
                    cellRealArray *f1 = static_cast<cellRealArray *>(Iter1->second);
                    setOutputField(nam, calcTotalPressure(
                                            "Ma_Ave", calcTimeAveVar(p()), calcTimeAveVar(gama()),
                                            calcTimeAveVar(v()), calcTimeAveVar(rho())) -
                                            *f1);
                } else if (nam == "Tt_Pulse") {
                    cellRealArray *f1 = static_cast<cellRealArray *>(Iter1->second);
                    if (flows_.mixtures().isSingular()) {
                        setOutputField(
                            nam, calcTotalTemperature("Ma_Ave", calcTimeAveVar(T()),
                                                      calcTimeAveVar(gama()), calcTimeAveVar(v()),
                                                      calcTimeAveVar(p()), calcTimeAveVar(rho())) -
                                     *f1);
                    } else {
                        calcTimeAveSpecies();
                        setOutputField("Tt_Pulse",
                                       calculateFieldVar::totalTemperature(
                                           flows_, calcTimeAveVar(p()), calcTimeAveVar(T()),
                                           calcTimeAveVar(v()), yiAve()) -
                                           *f1);
                    }
                } else {
                    if (Iter1->second->nElements() == 1) {
                        cellRealArray *f = static_cast<cellRealArray *>(Iter1->second);
                        setOutputField(nam, calcTimeAveVar(f->getTimeSumPtr()) - *f);
                    }
                    if (Iter1->second->nElements() == 3) {
                        cellVectorArray *f = static_cast<cellVectorArray *>(Iter1->second);
                        setOutputField(nam, calcTimeAveVar(f->getTimeSumPtr()) - *f);
                    }
                    break;
                }
            } else if (Iter1 == mesh().table().end() && Iter2 == p().tb().table().end()) {
                if (nam == "Ma_Pulse") {
                    setOutputField(nam,
                                   calculateFieldVar::calcMachNumber(
                                       calcTimeAveVar(v()), calcTimeAveVar(gama()),
                                       calcTimeAveVar(p()), calcTimeAveVar(rho())) -
                                       calculateFieldVar::calcMachNumber(v(), gama(), p(), rho()));
                }
                if (nam == "pt_Pulse") {
                    setOutputField(nam, calcTotalPressure(
                                            "Ma_Ave", calcTimeAveVar(p()), calcTimeAveVar(gama()),
                                            calcTimeAveVar(v()), calcTimeAveVar(rho())) -
                                            calcTotalPressure("Ma", p(), gama(), v(), rho()));
                }
                if (nam == "Tt_Pulse") {
                    if (flows_.mixtures().isSingular()) {
                        setOutputField(
                            nam, calcTotalTemperature("Ma_Ave", calcTimeAveVar(T()),
                                                      calcTimeAveVar(gama()), calcTimeAveVar(v()),
                                                      calcTimeAveVar(p()), calcTimeAveVar(rho())) -
                                     calcTotalTemperature("Ma", T(), gama(), v(), p(), rho()));
                    } else {
                        calcTimeAveSpecies();
                        setOutputField(nam, calculateFieldVar::totalTemperature(
                                                flows_, calcTimeAveVar(p()), calcTimeAveVar(T()),
                                                calcTimeAveVar(v()), yiAve()) -
                                                calculateFieldVar::totalTemperature(
                                                    flows_, p(), T(), v(), yi()));
                    }
                }
            }
        }
    }
}

const OpenHurricane::realArray
OpenHurricane::writePulseFieldVar::calcTotalPressure(const string &Ma, const realArray &P,
                                                     const realArray &gama, const vectorArray &v,
                                                     const realArray &rho) {
    auto iter = mesh().table().find(Ma);
    if (iter != mesh().table().end()) {
        cellRealArray *f = static_cast<cellRealArray *>(iter->second);
        return calculateFieldVar::calcTotalPressure(P, gama, *f, flows_);
    } else {
        return calculateFieldVar::calcTotalPressure(P, gama, v, rho, flows_);
    }
}

const OpenHurricane::realArray
OpenHurricane::writePulseFieldVar::calcTotalTemperature(const string &Ma, const realArray &T,
                                                        const realArray &gama, const vectorArray &v,
                                                        const realArray &P, const realArray &rho) {
    auto iter = mesh().table().find(Ma);
    if (iter != mesh().table().end()) {
        cellRealArray *f = static_cast<cellRealArray *>(iter->second);
        return calculateFieldVar::calcTotalTemperature(T, gama, *f, flows_);
    } else {
        return calculateFieldVar::calcTotalTemperature(T, gama, v, P, rho, flows_);
    }
}

void OpenHurricane::writePulseFieldVar::correct() const {
    sumField(mesh());
}

void OpenHurricane::writePulseFieldVar::reseting() const {
    writeFieldVar::reseting();
    reset();
}

void OpenHurricane::writePulseFieldVar::calcTimeAveSpecies() {
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

void OpenHurricane::writePulseFieldVar::setCurretTimeGap() {
    curretTimeGap_ = iter_.pTStep().totalTime() - startTime_;
}

OpenHurricane::fileName OpenHurricane::writePulseFieldVar::getFileName() const {
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

bool OpenHurricane::writePulseFieldVar::writeNow() const {
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
