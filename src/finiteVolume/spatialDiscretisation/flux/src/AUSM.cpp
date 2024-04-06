/*!
 * \file AUSM.cpp
 * \brief Main subroutines for AUSM scheme.
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

#include "AUSM.hpp"
#include <cmath>

namespace OpenHurricane {
    createClassNameStr(AUSM, "AUSM");
    registerObjFty(upwindScheme, AUSM, controller);
} // namespace OpenHurricane

void OpenHurricane::AUSM::calcFlux(const real rhol, const real rhor, const vector &VL,
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

    // The kinetic energy of left state
    //real ekl = 0.5 * (ul * ul + vl * vl + wl * wl);

    // The kinetic energy of right state
    //real ekr = 0.5 * (ur * ur + vr * vr + wr * wr);

    // Total enthalpy

    /*!\brief Total enthalpy of left state.*/
    real hl = el + pl / rhol;
    //real hl = hal + ekl;

    /*!\briefTotal enthalpy of right state.*/
    real hr = er + pr / rhor;
    //real hr = har + ekr;

    //// Total energy

    //	/*!\brief Total energy of left state.*/
    //	real el = hl - pl / rhol;

    //	/*!\brief Total energy of right state.*/
    //	real er = hr - pr / rhor;

    // Mach number

    /*!\brief Mach number at left state.*/
    real MaL = unl / cl;

    /*!\brief Mach number at right state.*/
    real MaR = unr / cr;

    /*!\brief The splitting Mach number at left state.*/
    real MaLP;
    if (MaL >= 1.0) {
        MaLP = MaL;
    } else if (fabs(MaL) < 1.0) {
        MaLP = 0.25 * sqr(MaL + 1.0);
    } else {
        MaLP = Zero;
    }

    /*!\brief The splitting Mach number at right state.*/
    real MaRM;
    if (MaR >= 1.0) {
        MaRM = Zero;
    } else if (fabs(MaR) < 1.0) {
        MaRM = -0.25 * sqr(MaR - 1.0);
    } else {
        MaRM = MaR;
    }

    real ManLR = MaLP + MaRM;

    if (flux.size() < 5) {
        flux.resize(5, Zero);
    }

    real plp;
    if (MaL >= 1.0) {
        plp = pl;
    } else if (fabs(MaL) < 1.0) {
        plp = 0.25 * pl * sqr(MaL + 1.0) * (2.0 - MaL);
    } else {
        plp = Zero;
    }

    real prm;
    if (MaR >= 1.0) {
        prm = Zero;
    } else if (fabs(MaR) < 1.0) {
        prm = 0.25 * pr * sqr(MaR - 1.0) * (2.0 + MaR);
    } else {
        prm = pr;
    }
    real plr = plp + prm;
    /*realArray pflux(5, Zero);
    pflux[1] = nx * plr;
    pflux[2] = ny * plr;
    pflux[3] = nz * plr;*/

    pFlux_[0] = 0.0;
    pFlux_[1] = nx * plr;
    pFlux_[2] = ny * plr;
    pFlux_[3] = nz * plr;
    pFlux_[4] = 0.0;

    if (ManLR >= 0.0) {
        real fmass = ManLR * rhol * cl;
        flux[0] = 1.0;
        flux[1] = ul;
        flux[2] = vl;
        flux[3] = wl;
        flux[4] = hl;
        flux *= fmass;
    } else {
        real fmass = ManLR * rhor * cr;
        flux[0] = 1.0;
        flux[1] = ur;
        flux[2] = vr;
        flux[3] = wr;
        flux[4] = hr;
        flux *= fmass;
    }

    //flux += pflux;
    flux += pFlux_;
    flux *= areaMag;
}

OpenHurricane::AUSM::AUSM(const controller &cont, const runtimeMesh &mesh, flowModel &flow)
    : upwindScheme(cont, mesh, flow), pFlux_(5, Zero) {}
