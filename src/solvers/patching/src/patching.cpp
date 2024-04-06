/*!
 * \file patching.cpp
 * \brief Main subroutines of class of patching.
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

#include "patching.hpp"

OpenHurricane::patching::patching(const runtimeMesh &mesh, const controller &cont)
    : mesh_(mesh), cont_(cont), regionName_(), typeName_(), varName_(), startStep_(0),
      stayStep_(0) {
    if (!cont_.found("region")) {
        errorAbortStr(("Cannot find option: \"region\" in " + cont_.name()));
    }
    regionName_ = cont_.findWord("region");

    if (!cont_.found("type")) {
        errorAbortStr(("Cannot find option: \"type\" in " + cont_.name()));
    }
    typeName_ = cont_.findWord("type");

    if (!cont_.found("patchVar")) {
        errorAbortStr(("Cannot find option: \"patchVar\" in" + cont_.name()));
    }

    std::string rnl = cont_.findText("patchVar");
    replaceAllMarks(rnl, "\n", " ");
    if (!rnl.empty()) {
        split(rnl, varName_, ",");
        for (integer i = 0; i < varName_.size(); ++i) {
            trim(varName_[i]);
        }
    }

    if (cont_.found("startStep")) {
        startStep_ = cont_.findType<integer>("startStep", integer(0));
    }
    if (cont_.found("stay")) {
        stayStep_ = cont_.findType<integer>("stay", integer(0));
    }
}

void OpenHurricane::patching::uniformPatch() const {
    const auto &varTable = mesh_.table();

    for (integer i = 0; i < varName_.size(); ++i) {
        auto iter = varTable.find(varName_[i]);
        if (iter != varTable.end()) {
            object *ob = iter->second;
            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                if (cont_.found(varName_[i])) {
                    real val = cont_.findType<real>(varName_[i], real(0.0));
                    mesh_.region(regionName_).patching(*f, val);
                }
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                if (cont_.found(varName_[i])) {
                    vector val = cont_.findType<vector>(varName_[i], vector(0.0));
                    mesh_.region(regionName_).patching(*f, val);
                }
            }
        }
    }
}

void OpenHurricane::patching::distributePatch() const {
    const auto &varTable = mesh_.table();

    for (integer i = 0; i < varName_.size(); ++i) {
        auto iter = varTable.find(varName_[i]);
        if (iter != varTable.end()) {
            object *ob = iter->second;
            if (ob->nElements() == 1) {
                cellRealArray *f = static_cast<cellRealArray *>(ob);
                if (cont_.found(varName_[i])) {
                    std::string formula = cont_.findText(varName_[i]);
                    replaceAllMarks(formula, "\n", " ");
                    mesh_.region(regionName_).distributing(*f, formula);
                }
            } else if (ob->nElements() == 3) {
                cellVectorArray *f = static_cast<cellVectorArray *>(ob);
                if (cont_.found(varName_[i])) {
                    std::string formula = cont_.findText(varName_[i]);
                    replaceAllMarks(formula, "\n", " ");
                    trim(formula);
                    mesh_.region(regionName_).distributing(*f, formula);
                }
            }
        }
    }
}
