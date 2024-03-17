/*!
 * \file monitor.cpp
 * \brief Main subroutines for monitor.
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

#include "monitor.hpp"
#include "controllerSwitch.hpp"

namespace OpenHurricane {
    createClassName(monitor);
    createObjFty(monitor, controller);
} // namespace OpenHurricane

void OpenHurricane::monitor::removeDuplicate(stringList &residualsNameList) const {
    if (residualsNameList.size() == 0) {
        return;
    }
    integer m = 0;
    stringList tmpList(residualsNameList.size());

    for (integer i = 0; i < residualsNameList.size(); ++i) {
        integer j = 0;
        for (; j < m; ++j) {
            if (residualsNameList[i] == tmpList[j]) {
                break;
            }
        }

        if (j == m) {
            tmpList[m] = residualsNameList[i];
            m++;
        }
    }
    tmpList.resize(m);
    residualsNameList.transfer(tmpList);
}

OpenHurricane::monitor::monitor(const iteration &iter, const runtimeMesh &mesh,
                                const controller &cont, const string &name)
    : iter_(iter), mesh_(mesh), updateStep_(cont.findOrDefault<integer>("updateStep", 1)),
      printToScreen_(true), writeToFile_(false), monitorName_(name), outFile_() {
    printToScreen_ = controllerSwitch(cont)("printToScreen", printToScreen_);

    writeToFile_ = controllerSwitch(cont)("writeToFile", writeToFile_);

    if (writeToFile_) {
        if (cont.found("fileName")) {
            string pw = cont.findWord("fileName");
            trim(pw);
            outFile_ = pw;
            if (!outFile_.isAbsolute()) {
                outFile_ = iter_.outputPath() / outFile_;
            }
        } else {
            outFile_ = iter_.configName().name(true) + "_" + monitorName_ + ".dat";
            outFile_ = iter_.outputPath() / outFile_;
        }
    }
}

hur_nodiscard OpenHurricane::uniquePtr<OpenHurricane::monitor>
OpenHurricane::monitor::creator(const iteration &iter, const runtimeMesh &mesh,
                                const controller &cont, const string &name) {
    string monitorType = cont.findWord("monitorType");
    defineInObjCreator(monitor, monitorType, controller, (iter, mesh, cont, name));
}