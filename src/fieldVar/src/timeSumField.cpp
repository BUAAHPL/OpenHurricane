/*!
 * \file timeSumField.cpp
 * \brief Main subroutines for timeSumField.
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
#include "timeSumField.hpp"

OpenHurricane::stringList OpenHurricane::timeSumField::timeSumVars_(0);
bool OpenHurricane::timeSumField::updating_(false);

void OpenHurricane::timeSumField::setTimeSumVarList(const flowModel &flows, const runtimeMesh &mesh,
                                                    const stringList &varList) const {
    for (integer i = 0; i < varList.size(); ++i) {
        if (mesh.foundOnlyObject(getPrefix(varList[i]))) {
            integer j = 0;
            for (j = 0; j < timeSumVars_.size(); ++j) {
                if (getPrefix(varList[i]) == timeSumVars_[j]) {
                    break;
                }
            }
            if (j == timeSumVars_.size()) {
                timeSumVars_.append(getPrefix(varList[i]));
            }
        } else if (getPrefix(varList[i]) == "Ma") {
            stringList varnameList = {"v", "gama", "p", "rho"};
            addTimeSumPriVar(varnameList);
        } else if (getPrefix(varList[i]) == "Pt") {
            stringList varnameList = {"p", "v", "gama", "p", "rho"};
            addTimeSumPriVar(varnameList);
        } else if (getPrefix(varList[i]) == "Tt") {
            if (flows.mixtures().isSingular()) {
                stringList varnameList = {"T", "v", "gama", "p", "rho"};
                addTimeSumPriVar(varnameList);
            } else {
                stringList varnameList;
                for (integer j = 0; j < flows.mixtures().species().size(); ++j) {
                    varnameList.append(flows.mixtures().species()[i].name());
                }
                addTimeSumPriVar(varnameList);
            }
        }
    }
}

void OpenHurricane::timeSumField::addTimeSumPriVar(const stringList &varnameList) const {
    for (integer i = 0; i < varnameList.size(); ++i) {
        integer j = 0;
        for (j = 0; j < timeSumVars_.size(); ++j) {
            if (varnameList[i] == timeSumVars_[j]) {
                break;
            }
        }
        if (j == timeSumVars_.size()) {
            timeSumVars_.append(varnameList[i]);
        }
    }
}

OpenHurricane::timeSumField::timeSumField() {}

OpenHurricane::timeSumField::~timeSumField() noexcept {}

void OpenHurricane::timeSumField::sumField(const runtimeMesh &mesh) const {
    const auto &iter = mesh.Iteration();

    if (iter.hasPhysicalTimeStep() && !updating_) {
        for (integer i = 0; i < timeSumVars_.size(); ++i) {
            if (mesh.foundOnlyObject(timeSumVars_[i])) {
                const auto &ob = mesh.findOnlyObject(timeSumVars_[i]);

                ob.calcTimeSumPtr(iter.pTStep().pTimeStep());
            }
#ifdef HUR_DEBUG
            else {
                LFatal("Can not find \"%s\" in mesh table", timeSumVars_[i].c_str());
            }
#endif // HUR_DEBUG
        }
        updating_ = true;
    }
}

bool OpenHurricane::timeSumField::reset() const {
    bool tmp = updating_;
    updating_ = false;
    return tmp;
}
