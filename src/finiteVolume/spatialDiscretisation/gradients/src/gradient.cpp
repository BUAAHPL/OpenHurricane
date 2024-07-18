/*!
 * \file gradient.cpp
 * \brief Main subroutines for gradient.
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

#include "gradient.hpp"
#include "cellMesh.hpp"
#include "faceMesh.hpp"
#include "transfer.hpp"

namespace OpenHurricane {
    createClassNameStr(gradient, "gradient");
    createObjFty(gradient, controller);
} // namespace OpenHurricane

OpenHurricane::uniquePtr<OpenHurricane::gradient>
OpenHurricane::gradient::creator(const controller &cont) {
    string gradientType = cont.findWord(gradient::className_);

    defineInObjCreator(gradient, static_cast<std::string>(gradientType), controller, (cont));
}

void OpenHurricane::gradient::updateBoundaryField(geometryArray<vector, cellMesh> &grad) {
    const auto &mesh = grad.mesh();

    const auto &fZ = mesh.faceZones();

    const auto &fL = mesh.faces();

    const auto &fS = mesh.faceArea();
    vectorTransfer myTransfer(mesh, grad, false, true);
    myTransfer.transferInit();
    for (integer fZI = 0; fZI < fZ.size(); ++fZI) {
        if (fZ[fZI].isInterior()) {
            continue;
        } else if (fZ[fZI].isSymmetric()) {
            for (integer fI = fZ[fZI].firstIndex(); fI <= fZ[fZI].lastIndex(); ++fI) {
                const auto &cl = fL[fI].leftCell();
                const auto &cr = fL[fI].rightCell();

                const vector n = fS[fI].normalized();

                grad[cr] = grad[cl] - real(2.0) * (n * grad[cl]) * n;
            }
        } else if (fZ[fZI].isCutFace()) {
            continue;
        } else if (fZ[fZI].isPeriodic() || fZ[fZI].isPeriodicShadow()) {
            continue;
        } else {
            for (integer fI = fZ[fZI].firstIndex(); fI <= fZ[fZI].lastIndex(); ++fI) {
                const auto &cl = fL[fI].leftCell();
                const auto &cr = fL[fI].rightCell();
                grad[cr] = grad[cl];
            }
        }
    }
    myTransfer.transferring();
}

void OpenHurricane::gradient::updateBoundaryField(geometryArray<tensor, cellMesh> &grad) {
    const auto &mesh = grad.mesh();

    const auto &fZ = mesh.faceZones();

    const auto &fL = mesh.faces();

    const auto &fS = mesh.faceArea();
    tensorTransfer myTransfer(mesh, grad, false, true);
    myTransfer.transferInit();

    tensor nfnf;
    tensor HousHolder;

    for (integer fZI = 0; fZI < fZ.size(); ++fZI) {
        if (fZ[fZI].isInterior()) {
            continue;
        } else if (fZ[fZI].isSymmetric()) {
            //for (integer fI = fZ[fZI].firstIndex(); fI <= fZ[fZI].lastIndex(); ++fI) {
            //    const auto &cl = fL[fI].leftCell();
            //    const auto &cr = fL[fI].rightCell();

            //    const vector n = fS[fI].normalized();

            //    grad[cr] = grad[cl] - (real(2.0) * (grad[cl] * n) & n);
            //}

            for (integer fI = fZ[fZI].firstIndex(); fI <= fZ[fZI].lastIndex(); ++fI) {
                const auto &cl = fL[fI].leftCell();
                const auto &cr = fL[fI].rightCell();

                const vector n = fS[fI].normalized();

                /*grad[cr] = grad[cl] - (real(2.0) * (grad[cl] * n) & n);*/

                for (integer i = 0; i < 3; i++) {
                    for (integer j = 0; j < 3; j++) {
                        nfnf(i,j) = n[i] * n[j];
                    }
                }

                HousHolder = I - real(2) * nfnf;

                grad[cr] = HousHolder * grad[cl] * HousHolder;
            }
        } else if (fZ[fZI].isCutFace()) {
            continue;
        } else if (fZ[fZI].isPeriodic() || fZ[fZI].isPeriodicShadow()) {
            continue;
        } else {
            for (integer fI = fZ[fZI].firstIndex(); fI <= fZ[fZI].lastIndex(); ++fI) {
                const auto &cl = fL[fI].leftCell();
                const auto &cr = fL[fI].rightCell();
                grad[cr] = grad[cl];
            }
        }
    }

    myTransfer.transferring();
}
