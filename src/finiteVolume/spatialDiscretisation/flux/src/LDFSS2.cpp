/*!
 * \file LDFSS2.cpp
 * \brief Main subroutines for LDFSS2 scheme.
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

#include "LDFSS2.hpp"
#include <cmath>

const OpenHurricane::real OpenHurricane::LDFSS2::delta_(1.0);

namespace OpenHurricane {
    createClassNameStr(LDFSS2, "LDFSS2");
    registerObjFty(upwindScheme, LDFSS2, controller);
} // namespace OpenHurricane

void OpenHurricane::LDFSS2::calcFlux(const real rhol, const real rhor, const vector &VL,
                                     const vector &VR, const real pl, const real pr, const real gl,
                                     const real gr, const real el, const real er, const real cl,
                                     const real cr, const vector &faceArea, const real blend,
                                     realArray &flux) const {
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

    //// Total energy

    //	// Total energy of left state
    //	real el = hl - pl / rhol;

    //	// Total energy of right state
    //	real er = hr - pr / rhor;

    // Mach number

    // The acoustic speed at the interface
    const real cLR = 0.5 * (cl + cr);

    /*!\brief Mach number at left state.*/
    const real MaL = unl / cLR;

    /*!\brief Mach number at right state.*/
    const real MaR = unr / cLR;

    /*!\brief The splitting Mach number at left state.*/
    const real MaLP = 0.25 * sqr(MaL + real(1.0));

    /*!\brief The splitting Mach number at right state.*/
    const real MaRM = -0.25 * sqr(MaR - real(1.0));

    const real betaL = -max(Zero, real(1 - int(mag(MaL))));
    const real betaR = -max(Zero, real(1 - int(mag(MaR))));

    const real alphaLP = 0.5 * (1.0 + sign(MaL));
    const real alphaRM = 0.5 * (1.0 - sign(MaR));

    const real C_VLP = alphaLP * (1.0 + betaL) * MaL - betaL * MaLP;
    const real C_VRM = alphaRM * (1.0 + betaR) * MaR - betaR * MaRM;

    const real MaLR =
        0.25 * betaL * betaR * sqr(sqrt(real(0.5) * (sqr(MaL) + sqr(MaR))) - real(1.0));

    const real MaLRP = MaLR * (2.0 * pr / (pl + pr) - delta_ * mag(pl - pr) / pl);
    const real MaLRM = MaLR * (2.0 * pl / (pl + pr) - delta_ * mag(pl - pr) / pr);

    const real ManLRP = C_VLP - MaLRP;
    const real ManLRM = C_VRM + MaLRM;

    const real DLP = alphaLP * (1.0 + betaL) - betaL * P5P(MaL, real(0.0));
    const real DRM = alphaRM * (1.0 + betaR) - betaR * P5M(MaR, real(0.0));
    const real plr = DLP * pl + DRM * pr;

    if (flux.size() < 5) {
        flux.resize(5, Zero);
    }

    pFlux_[0] = 0.0;
    pFlux_[1] = nx * plr;
    pFlux_[2] = ny * plr;
    pFlux_[3] = nz * plr;
    pFlux_[4] = 0.0;

    const real fmassLP = rhol * cLR * ManLRP;
    const real fmassRM = rhor * cLR * ManLRM;

    massFluxL_[0] = 1.0;
    massFluxL_[1] = ul;
    massFluxL_[2] = vl;
    massFluxL_[3] = wl;
    massFluxL_[4] = hl;
    massFluxL_ *= fmassLP;

    massFluxR_[0] = 1.0;
    massFluxR_[1] = ur;
    massFluxR_[2] = vr;
    massFluxR_[3] = wr;
    massFluxR_[4] = hr;
    massFluxR_ *= fmassRM;

    //flux += pflux;
    flux = massFluxL_;
    flux += massFluxR_;
    flux += pFlux_;
    flux *= areaMag;
}

OpenHurricane::LDFSS2::LDFSS2(const controller &cont, const runtimeMesh &mesh, flowModel &flow)
    : upwindScheme(cont, mesh, flow), pFlux_(5), massFluxL_(5), massFluxR_(5) {}
