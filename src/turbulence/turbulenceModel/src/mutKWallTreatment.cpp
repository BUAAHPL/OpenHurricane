/*!
 * \file mutKWallTreatment.cpp
 * \brief Main subroutines of mut near-wall treatment.
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
#include "mutKWallTreatment.hpp"

namespace OpenHurricane {
    createClassNameStr(mutKWallTreatment, "mutKWallTreatment");
    registerObjFty(realBoundary, mutKWallTreatment, controller);
} // namespace OpenHurricane

OpenHurricane::real OpenHurricane::mutKWallTreatment::yPlusLam(const real kappa, const real E) {
    real yl = 11.0;
    for (uinteger i = 0; i < 10; ++i) {
        yl = log(max(E * yl, real(1))) / kappa;
    }
    return yl;
}

OpenHurricane::mutKWallTreatment::mutKWallTreatment(const faceZone &fZ,
                                                    geometryArray<real, cellMesh> &gf,
                                                    const controller &cont)
    : Base(fZ, gf, cont), beta_(3.0 / 40.0), Cmiu_(0.090), kappa_(0.41), E_(9.8), yplusLam_(11.0) {
    yplusLam_ = yPlusLam(kappa_, E_);
}

OpenHurricane::mutKWallTreatment::mutKWallTreatment(const mutKWallTreatment &bB)
    : Base(bB), beta_(bB.beta_), Cmiu_(bB.Cmiu_), kappa_(bB.kappa_), E_(bB.E_),
      yplusLam_(bB.yplusLam_) {}

void OpenHurricane::mutKWallTreatment::updateBoundary() {
    const registerTable &tb = this->varField().tb();
    const auto &y = tb.findObjectRef<cellRealArray>("wallDist");
    const auto &rho = tb.findObjectRef<cellRealArray>("rho");
    const auto &v = tb.findObjectRef<cellVectorArray>("v");
    const auto &kt = tb.findObjectRef<cellRealArray>("kt");
    const auto &wt = tb.findObjectRef<cellRealArray>("wt");
    const auto &mu = tb.findObjectRef<cellRealArray>("mu");

    const auto &faces = mesh().faces();
    const auto &fA = mesh().faceArea();
    auto &mut = varArray_;
    const real Cmu025 = pow(Cmiu_, real(0.25));
    for (integer facei = boundaryZone_.firstIndex(); facei <= boundaryZone_.lastIndex(); ++facei) {
        const integer cl = faces[facei].leftCell();
        const integer cr = faces[facei].rightCell();

        real rhow = (rho[cl] + rho[cr]) * real(0.5);
        real muw = (mu[cl] + mu[cr]) * real(0.50);
        const real yPlus = Cmu025 * y[cl] * sqrt(kt[cl]) / (muw / rhow);

        real mutw = 0;
        if (yPlus > yplusLam_) {
            // log - region
            mutw = muw * (yPlus * kappa_ / log(E_ * yPlus) - 1.0);
        } else {
            // laminar sub-layer
        }

        mut[cr] = 2.0 * mutw - mut[cl];
    }
}
