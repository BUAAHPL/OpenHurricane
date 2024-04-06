/*!
 * \file loadBalancingWeights.cpp
 * \brief Main subroutines for decomposing mesh.
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

#include "loadBalancingWeights.hpp"

OpenHurricane::loadBalancingWeights::loadBalancingWeights(const originMeshRead &originMeshes,
                                                          const controller &cont)
    : originMesh_(originMeshes), weightType_(weightTypes::unweight), faceAdditionalWeights_(2) {
    if (cont.found("weightType")) {
        const auto ww = cont.findWord("weightType");

        if (ww == "unweight") {
            weightType_ = weightTypes::unweight;
            Pout << "    Info: Load balancing weights using unweight..." << std::endl;
        } else if (ww == "facesPerCell") {
            weightType_ = weightTypes::facesPerCell;
            if (cont.found("facesPerCell")) {
                const auto &fpcCont = cont.subController("facesPerCell");
                faceAdditionalWeights_ =
                    fpcCont.findOrDefault<integer>("additionalWeights", faceAdditionalWeights_);
            }
            Pout << "    Info: Load balancing weights using faces-per-cell "
                    "with an additional weight: "
                 << faceAdditionalWeights_ << "..." << std::endl;
        } else if (ww == "givenInMeshFile") {
            weightType_ = weightTypes::givenInMeshFile;
            Pout << "    Info: Load balancing weights using weight given in "
                    "mesh file..."
                 << std::endl;
        } else {
            errorAbortStr(("Unknown weight type: " + ww + " in " + cont.name()));
        }
    }
}

hur_nodiscard OpenHurricane::integerArray OpenHurricane::loadBalancingWeights::weights() const {
    integerArray wgt;
    if (isUnWeight()) {
        wgt.resize(originMesh_.cells().size(), 1);
        return wgt;
    } else if (isFacesPerCell()) {
        wgt.resize(originMesh_.cells().size(), faceAdditionalWeights_);
        for (integer i = 0; i < originMesh_.cells().size(); ++i) {
            wgt[i] += originMesh_.cells()[i].faceSize();
        }
    } else if (isGivenByMeshFile()) {
        wgt.resize(originMesh_.cells().size(), 1);
        for (integer i = 0; i < originMesh_.cellLoadWeights().size(); ++i) {
            wgt[i] = originMesh_.cellLoadWeights()[i];
        }
    }
    return wgt;
}