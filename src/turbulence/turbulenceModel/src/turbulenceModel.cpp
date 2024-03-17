/*!
 * \file turbulenceModel.cpp
 * \brief Main subroutines for the turbulence model.
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

#include "turbulenceModel.hpp"

 namespace OpenHurricane{
	createClassNameStr(turbulenceModel,"turbulenceModel");
}

namespace OpenHurricane {
     createObjFty(turbulenceModel,controller);
    void turbulenceModel::solving(const realArray &dt) {}

    bool turbulenceModel::loop() {
        if (isCoupled()) {
            return false;
        }
        bool endding = cStep_ >= maxStep_;
        if (endding) {
            cStep_ = 0;
        } else {
            cStep_++;
        }
        return !endding;
    }

    faceSymmTensorArray turbulenceModel::tauEff(const faceRealArray &rhof,
                                                const faceRealArray &mulf,
                                                const faceRealArray &mutf,
                                                const faceTensorArray &deltafV) const {
        return faceSymmTensorArray::nullObject();
    }

    symmTensorArray turbulenceModel::ReynoldsStressTensor() const {
        LFatal("Only available for RANS models");
        return symmTensorArray();
    }
} // namespace OpenHurricane

OpenHurricane::turbulenceModel::turbulenceModel(const controller &cont, flowModel &flowMod)
    : yPtr_(nullptr), flowM_(flowMod), nEq_(0), solveType_(coupled), cStep_(0), maxStep_(1),
      mul_(flowMod.mul()), mut_(flowMod.mut()) {
    if (cont.found("solverType")) {
        string stw = cont.findWord("solverType");
        trim(stw);
        stringToUpperCase(stw);

        if (stw == "COUPLED") {
            solveType_ = coupled;
            maxStep_ = 0;
        } else if (stw == "SPLITTING") {
            solveType_ = splitting;
            //maxStep_ = cont.findOrDefault<integer>("maxStep", 1);
            maxStep_ = 1;
        } else {
            LFatal("Unknown solver type for turbulence");
        }
    }

}

OpenHurricane::uniquePtr<OpenHurricane::turbulenceModel>
OpenHurricane::turbulenceModel::creator(const controller &cont, flowModel &flowMod) {
    string turbulenceType = cont.findWord(turbulenceModel::className_);

    Pout << "    Setting turbulence model: " << turbulenceType << std::endl;
	defineInObjCreator(turbulenceModel,static_cast<std::string>(turbulenceType),controller,(cont, flowMod));
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::turbulenceModel::up(const realArray &muw, const integer zoneId) {
    const auto &fz = mesh().faceZones()[zoneId];
    const auto &fA = mesh().faceArea();
    const auto &faces = mesh().faces();

    realArray upw(fz.size());

    integer id = 0;
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const integer &cl = faces[facei].leftCell();

        real uu = v()[cl] * fA[facei].normalized();

        real ut = sqrt(v()[cl].magSqr() - sqr(uu));
        real tauw = muw[id] * ut / (wallDist()[cl]);
        real utau = sqrt(tauw / rho()[cl]);
        upw[id++] = ut / utau;
    }

    return upw;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::turbulenceModel::up(const integer zoneId) {
    return up(mul(zoneId), zoneId);
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::turbulenceModel::yPlus(const integer zoneId) {
    const auto &fz = mesh().faceZones()[zoneId];
    const auto &fA = mesh().faceArea();
    const auto &faces = mesh().faces();

    realArray ypw(fz.size());
    realArray muw = mul(zoneId);
    integer id = 0;
    for (integer facei = fz.firstIndex(); facei <= fz.lastIndex(); ++facei) {
        const auto &cl = faces[facei].leftCell();
        const auto &cr = faces[facei].rightCell();

        real uu = v()[cl] * fA[facei].normalized();

        real rhow = real(0.5) * (rho()[cl] + rho()[cr]);
        real ut = sqrt(v()[cl].magSqr() - sqr(uu));
        real tauw = muw[id] * ut / (wallDist()[cl]);
        real utau = sqrt(tauw / rhow);
        ypw[id] = wallDist()[cl] * utau / (muw[id] / rhow);
        id++;
    }
    return ypw;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::turbulenceModel::KolmogorovLengthScale() const {
    realArray yitaK(mesh().nTotalCells(), Zero);
    const auto nu = flowM_.nu();
    const auto ep = epsilon();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        yitaK[n] = pow025(pow3(nu[n]) / max(ep[n], veryTiny));
    }

    return yitaK;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::turbulenceModel::KolmogorovTimeScale() const {
    realArray taoK(mesh().nTotalCells(), Zero);
    const auto nu = flowM_.nu();
    const auto ep = epsilon();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        taoK[n] = sqrt(nu[n] / max(ep[n], veryTiny));
    }

    return taoK;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::turbulenceModel::KolmogorovVelocityScale() const {
    realArray uK(mesh().nTotalCells(), Zero);
    const auto nu = flowM_.nu();
    const auto ep = epsilon();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        uK[n] = pow025(nu[n] * ep[n]);
    }

    return uK;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::turbulenceModel::integralLengthScale() const {
    realArray yita0(mesh().nTotalCells(), Zero);
    const auto kk = k();
    const auto ep = epsilon();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        yita0[n] = kk[n] * sqrt(kk[n]) / max(ep[n], veryTiny);
    }

    return yita0;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::turbulenceModel::integralTimeScale() const {
    realArray tao0(mesh().nTotalCells(), Zero);
    const auto kk = k();
    const auto ep = epsilon();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        tao0[n] = kk[n] / max(ep[n], veryTiny);
    }

    return tao0;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::turbulenceModel::integralVelocityScale() const {
    realArray u0(mesh().nTotalCells(), Zero);
    const auto kk = k();
    for (integer n = 0; n < mesh().nCells(); ++n) {
        u0[n] = sqrt(kk[n]);
    }

    return u0;
}

