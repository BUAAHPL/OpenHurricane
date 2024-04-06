/*!
 * \file AUSMPWPlus.cpp
 * \brief Main subroutines for AUSMPW+ scheme.
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

#include "AUSMPWPlus.hpp"
#include <cmath>

namespace OpenHurricane {
    createClassNameStr(AUSMPWPlus, "AUSMPWPlus");
    registerObjFty(upwindScheme, AUSMPWPlus, controller);
} // namespace OpenHurricane

void OpenHurricane::AUSMPWPlus::calcFlux(const real rhol, const real rhor, const vector &VL,
                                         const vector &VR, const real pl, const real pr,
                                         const real gl, const real gr, const real el, const real er,
                                         const real cl, const real cr, const vector &faceArea,
                                         const real blend, realArray &flux) const {
    // Faace area (magnitude)
    const real areaMag = faceArea.magnitude();

    // Unit normal vector components
    const real nx = faceArea.x() / areaMag;
    const real ny = faceArea.y() / areaMag;
    const real nz = faceArea.z() / areaMag;

    const real ul = VL[0];
    const real ur = VR[0];
    const real vl = VL[1];
    const real vr = VR[1];
    const real wl = VL[2];
    const real wr = VR[2];

    // Contravariant velocity

    /*!\brief Left state contravariant velocity.*/
    const real unl = ul * nx + vl * ny + wl * nz;

    /*!\brief Right state contravariant velocity.*/
    const real unr = ur * nx + vr * ny + wr * nz;

    // kinetic energy

    /*!\brief The kinetic energy of left state.*/
    //real ekl = 0.5 * (ul * ul + vl * vl + wl * wl);

    /*!\brief The kinetic energy of right state.*/
    //real ekr = 0.5 * (ur * ur + vr * vr + wr * wr);

    // Total enthalpy

    /*!\brief Total enthalpy of left state.*/
    real hl = el + pl / rhol;
    //real hl = hal + ekl;

    /*!\brief Total enthalpy of right state.*/
    real hr = er + pr / rhor;
    //real hr = har + ekr;

    real cLR = OpenHurricane::min(cl, cr);
    if (0.5 * (unl + unr) >= 0.0) {
        cLR = cLR * cLR / max(unl, cLR);
    } else {
        cLR = cLR * cLR / max(unr, cLR);
    }

    real rhoLR = 0.5 * (rhol + rhor);

    // Mach number

    /*!\brief Mach number at left state.*/
    real MaL = unl / cLR;

    /*!\brief Mach number at right state.*/
    real MaR = unr / cLR;

    real MaLP;
    real psiLP;
    if (std::abs(MaL) < 1.0) {
        MaLP = 0.25 * OpenHurricane::sqr(MaL + 1.0);
        psiLP = MaLP * (2.0 - MaL);
    } else {
        MaLP = 0.5 * (MaL + std::abs(MaL));
        psiLP = MaLP / MaL;
    }

    real MaRM;
    real psiRM;
    if (std::abs(MaR) < 1.0) {
        MaRM = -0.25 * OpenHurricane::sqr(MaR - 1.0);
        psiRM = -MaRM * (2.0 + MaR);
    } else {
        MaRM = 0.5 * (MaR - std::abs(MaR));
        psiRM = MaRM / MaR;
    }

    real MaLR = MaLP + MaRM;
    real PLR = psiLP * pl + psiRM * pr;

    real omega = min(pr / pl, pl / pr);
    omega = 1.0 - std::pow(omega, 3.0);
    real fL = Zero;
    if (std::abs(MaL) < 1.0) {
        fL = pl / PLR - 1.0;
    }
    real fR = Zero;
    if (std::abs(MaR) < 1.0) {
        fR = pr / PLR - 1.0;
    }

    real MLP;
    real MRM;

    if (MaLR >= 0.0) {
        MLP = MaLP + MaRM * ((1.0 - omega) * (1.0 + fR) - fL);
        MRM = MaRM * omega * (1.0 + fR);
    } else {
        MLP = MaLP * omega * (1.0 + fL);
        MRM = MaRM + MaLP * ((1.0 - omega) * (1.0 + fL) - fR);
    }

    /*!\brief mass flux.*/
    real massLP = MLP * cLR * rhol;
    real massRM = MRM * cLR * rhor;

    if (flux.size() < 5) {
        flux.resize(5, Zero);
    }

    flux[0] = 1.0;
    flux[1] = ul;
    flux[2] = vl;
    flux[3] = wl;
    flux[4] = hl;
    flux *= massLP;

    /*realArray fluxR(5, Zero);
    fluxR[0] = 1.0;
    fluxR[1] = ur;
    fluxR[2] = vr;
    fluxR[3] = wr;
    fluxR[4] = hr;
    fluxR *= massRM;*/

    fluxR_[0] = 1.0;
    fluxR_[1] = ur;
    fluxR_[2] = vr;
    fluxR_[3] = wr;
    fluxR_[4] = hr;
    fluxR_ *= massRM;

    //flux += fluxR;
    flux += fluxR_;
    /*!\brief Pressure flux.*/
    /*realArray fp(5, Zero);
    fp[1] = nx;
    fp[2] = ny;
    fp[3] = nz;
    fp = PLR * fp;
    flux += fp;*/

    fp_[0] = 0.0;
    fp_[1] = nx;
    fp_[2] = ny;
    fp_[3] = nz;
    fp_[4] = 0.0;
    fp_ = PLR * fp_;
    flux += fp_;

    flux *= areaMag;
}

void OpenHurricane::AUSMPWPlus::calcAS(const vector &VL, const vector &VR, const real hsl,
                                       const real hsr, const real gamml, const real gammr,
                                       const vector &faceArea, real &al, real &ar) const {
    // Faace area (magnitude)
    const real areaMag = faceArea.magnitude();

    // Unit normal vector components
    const real nx = faceArea.x() / areaMag;
    const real ny = faceArea.y() / areaMag;
    const real nz = faceArea.z() / areaMag;

    const real ul = VL[0];
    const real ur = VR[0];
    const real vl = VL[1];
    const real vr = VR[1];
    const real wl = VL[2];
    const real wr = VR[2];

    // Contravariant velocity

    /*!\brief Left state contravariant velocity.*/
    const real unl = ul * nx + vl * ny + wl * nz;

    /*!\brief Right state contravariant velocity.*/
    const real unr = ur * nx + vr * ny + wr * nz;

    /*real ekl = 0.5 * (sqr(ul) + sqr(vl) + sqr(wl));
    real ekr = 0.5 * (sqr(ur) + sqr(vr) + sqr(wr));*/

    real hnormall = hsl + 0.5 * unl * unl;
    real hnormalr = hsr + 0.5 * unr * unr;

    real gammaLRp1 = 0.5 * (gamml + gammr) + 1.0;
    real gammaLRm1 = 0.5 * (gamml + gammr) - 1.0;
    al = sqrt(2.0 * hnormall * gammaLRm1 / gammaLRp1);
    ar = sqrt(2.0 * hnormalr * gammaLRm1 / gammaLRp1);
}

OpenHurricane::AUSMPWPlus::AUSMPWPlus(const controller &cont, const runtimeMesh &mesh,
                                      flowModel &flow)
    : upwindScheme(cont, mesh, flow), fp_(5, Zero), fluxR_(5, Zero) {}
