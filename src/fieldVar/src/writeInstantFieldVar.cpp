/*!
 * \file writeInstantFieldVar.cpp
 * \brief Main subroutines for writing instant field variables.
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

#include "writeInstantFieldVar.hpp"
#include "calculateFieldVar.hpp"

namespace OpenHurricane {
    createClassNameStr(writeInstantFieldVar, "writeInstantFieldVar");
    registerObjFty(writeFieldVar, writeInstantFieldVar, controller);
} // namespace OpenHurricane

OpenHurricane::writeInstantFieldVar::writeInstantFieldVar(
    const flowModel &flows, const iteration &iter, const controller &cont, const string &writeId,
    std::map<std::string, object *> &outFieldVarMap)
    : writeFieldVar(flows, iter, cont, writeId, outFieldVarMap) {
    setoutVarMap(writeFieldVarList_);
}

void OpenHurricane::writeInstantFieldVar::updating() {
    if (!writeNow()) {
        return;
    }

    calcStagnationParameters();
    calcViscousRatio();
    calcVorticity();
    calcDeltaVorticity();
    calcOtherVorticity();
    calcMoleSpecies();
}

OpenHurricane::fileName OpenHurricane::writeInstantFieldVar::getFileName() const {
    fileName outN = iter_.outputName();
    auto pathOut = outN.parentPath();
    auto fname = outN.name(true);
    string fext;
    fext = ".plt";

    fname += "-";
    fname += "Instant";
    fname += "-";
    fname += toString(iter_.totalStep());
    outN = fname + fext;
    outN = pathOut / outN;

    return outN;
}

bool OpenHurricane::writeInstantFieldVar::writeNow() const {
    if (!isWritting_ || writeFieldVarList_.size() == 0) {
        return false;
    } else {
        if (iter_.cStep() % iter_.writeOutputStep() == 0 ||
            (iter_.end() && !argParse::noWriteAtEnd())) {
            return true;
        }
    }
    return false;
}

void OpenHurricane::writeInstantFieldVar::calcStagnationParameters() const {
    const auto iter = sets_.find("Ma");
    if (iter != sets_.end()) {
        setOutputField(*iter, calculateFieldVar::calcMachNumber(v(), gama(), p(), rho()));
    }
    const auto iter2 = sets_.find("pt");
    if (iter2 != sets_.end()) {
        setOutputField(*iter2, calculateFieldVar::calcTotalPressure(flows_));
    }

    const auto iterPtCp = sets_.find("pt_constCp");
    if (iterPtCp != sets_.end()) {
        setOutputField(*iterPtCp,
                       calculateFieldVar::calcTotalPressure(p(), gama(), v(), rho(), flows_));
    }

    const auto iter3 = sets_.find("Tt");
    if (iter3 != sets_.end()) {
        if (flows_.mixtures().isSingular()) {
            calcTotalTemperature("Ma", "Tt", T(), gama(), v(), p(), rho());
        } else {
            setOutputField(*iter3,
                           calculateFieldVar::totalTemperature(flows_, p(), T(), v(), yi()));
        }
    }
    const auto iterTtCp = sets_.find("Tt_constCp");
    if (iterTtCp != sets_.end()) {
        setOutputField(*iterTtCp, calculateFieldVar::calcTotalTemperature(T(), gama(), v(), p(),
                                                                          rho(), flows_));
    }

    const auto iterHat = sets_.find("hat");
    if (iterHat != sets_.end()) {
        setOutputField(*iterHat, calculateFieldVar::calcTotalEnthalpy(flows_));
    }

    const auto iterCV = sets_.find("cellVolume");
    if (iterCV != sets_.end()) {
        setOutputField(*iterCV, calculateFieldVar::cellVolume(flows_));
    }
}

void OpenHurricane::writeInstantFieldVar::calcViscousRatio() const {
    const auto iter = sets_.find("viscousRatio");
    if (iter != sets_.end()) {
        setOutputField(*iter, calculateFieldVar::calcViscousRatio(flows_));
    }
}

void OpenHurricane::writeInstantFieldVar::calcVorticity() const {
    const auto iter = sets_.find("vorticity");
    if (iter != sets_.end()) {
        setOutputField(*iter, calculateFieldVar::calcVorticity(flows_));
    }

    const auto iterQ = sets_.find("QCriterion");
    if (iterQ != sets_.end()) {
        setOutputField(*iterQ, calculateFieldVar::calcQCriterion(flows_));
    }

    const auto iterO = sets_.find("OmegaCriterion");
    if (iterO != sets_.end()) {
        setOutputField(*iterO, calculateFieldVar::calcOmegaCriterion(flows_));
    }

    const auto iterDel = sets_.find("DeltaCriterion");
    if (iterDel != sets_.end()) {
        setOutputField(*iterDel, calculateFieldVar::calcDeltaCriterion(flows_));
    }
}

void OpenHurricane::writeInstantFieldVar::calcDeltaVorticity() const {
    const auto iter = sets_.find("DeltaVorticity");
    if (iter != sets_.end()) {
        realArray P(-tr(v().grad()));
        realArray R(det(v().grad()));
        realArray Q(0.5 * (P * P - tr(v().grad() * v().grad())));
        setOutputField(*iter, Q * Q * Q / 27 + R * R / 4);
    }
}

void OpenHurricane::writeInstantFieldVar::calcOtherVorticity() const {
    const auto iter1 = sets_.find("QVorticity");
    if (iter1 != sets_.end()) {
        tensorArray Qtensor(skew(v().grad()) * skew(v().grad()) +
                            symm(v().grad()) * symm(v().grad()));
        if (iter1 != sets_.end()) {
            setOutputField(*iter1, -0.5 * tr(Qtensor));
        }
    }
}

void OpenHurricane::writeInstantFieldVar::calcMoleSpecies() const {
    if (flows_.mixtures().isSingular()) {
        return;
    }
    bool isOutXi = false;
    for (integer i = 0; i < flows_.mixtures().species().size(); ++i) {
        auto iter = sets_.find(flows_.mixtures().species()[i].name() + "-mole");
        if (iter != sets_.end()) {
            isOutXi = true;
            break;
        }
    }
    if (!isOutXi) {
        return;
    }
    PtrList<cellRealArray> tmpXi;
    for (integer i = 0; i < flows_.mixtures().species().size(); ++i) {
        tmpXi.append(new cellRealArray(
            object(flows_.mixtures().species()[i].name() + "-mole", mesh(), object::WRITE_OUTPUT),
            mesh()));
    }
    flows_.mixtures().species().Yi2Xi(yi(), tmpXi);
    for (integer i = 0; i < flows_.mixtures().species().size(); ++i) {
        auto iter = outFieldVarMap_.find(flows_.mixtures().species()[i].name() + "-mole");
        if (iter != outFieldVarMap_.end()) {
            iter->second = tmpXi.set(i, nullptr);
        }
    }
}