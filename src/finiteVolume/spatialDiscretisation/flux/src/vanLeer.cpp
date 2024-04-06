/*!
 * \file HLL.cpp
 * \brief Main subroutines for van Leer scheme.
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

#include "vanLeer.hpp"
#include <cmath>

namespace OpenHurricane {
    createClassNameStr(vanLeer, "vanLeer");
    registerObjFty(upwindScheme, vanLeer, controller);
} // namespace OpenHurricane

void OpenHurricane::vanLeer::calcFlux(const real rhol, const real rhor, const vector &VL,
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
    real ekl = 0.5 * (ul * ul + vl * vl + wl * wl);

    /*!\brief The kinetic energy of right state.*/
    real ekr = 0.5 * (ur * ur + vr * vr + wr * wr);

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

    /*!\brief Mach number at left state.*/
    real MaL = unl / cl;

    /*!\brief Mach number at right state.*/
    real MaR = unr / cr;

    /*!\brief The splitting Mach number at left state.*/
    real MaLP;
    if (MaL >= 1.0) {
        MaLP = MaL;
    } else if (std::abs(MaL) < 1.0) {
        MaLP = 0.25 * OpenHurricane::sqr(MaL + 1.0);
    } else {
        MaLP = Zero;
    }

    /*!\brief The splitting Mach number at right state.*/
    real MaRM;
    if (MaR >= 1.0) {
        MaRM = Zero;
    } else if (std::abs(MaR) < 1.0) {
        MaRM = -0.25 * OpenHurricane::sqr(MaR - 1.0);
    } else {
        MaRM = MaR;
    }

    real ManLR = MaLP + MaRM;

    if (flux.size() < 5) {
        flux.resize(5, Zero);
    }

    if (ManLR >= 1.0) {
        real rhou = rhol * unl;
        flux[0] = areaMag * rhou;
        flux[1] = areaMag * (rhou * ul + pl * nx);
        flux[2] = areaMag * (rhou * vl + pl * ny);
        flux[3] = areaMag * (rhou * wl + pl * nz);
        flux[4] = areaMag * rhou * hl;
    } else if (ManLR <= -1.0) {
        real rhou = rhor * unr;
        flux[0] = areaMag * rhou;
        flux[1] = areaMag * (rhou * ur + pr * nx);
        flux[2] = areaMag * (rhou * vr + pr * ny);
        flux[3] = areaMag * (rhou * wr + pr * nz);
        flux[4] = areaMag * rhou * hr;
    } else // |ManLR| < 1
    {
        real fmass = 0.25 * rhol * cl * OpenHurricane::sqr(MaL + 1.0);
        real mid = (2.0 * cl - unl) / gl;
        real fenergy = 0.5 * OpenHurricane::sqr((gl - 1.0) * unl + 2.0 * cl) / (gl * gl - 1.0) +
                       ekl - 0.5 * unl * unl;
        flux[0] = fmass;
        flux[1] = fmass * (nx * mid + ul);
        flux[2] = fmass * (ny * mid + vl);
        flux[3] = fmass * (nz * mid + wl);
        flux[4] = fmass * fenergy;

        fmass = -0.25 * rhor * cr * OpenHurricane::sqr(MaR - 1.0);
        mid = (-2.0 * cr - unr) / gr;
        fenergy = 0.5 * OpenHurricane::sqr((gr - 1.0) * unl - 2.0 * cr) / (gr * gr - 1.0) + ekr -
                  0.5 * unr * unr;
        flux[0] = flux[0] + fmass;
        flux[1] = flux[1] + fmass * (nx * mid + ur);
        flux[2] = flux[2] + fmass * (ny * mid + vr);
        flux[3] = flux[3] + fmass * (nz * mid + wr);
        flux[4] = flux[4] + fmass * fenergy;

        flux *= areaMag;
    }
}

OpenHurricane::vanLeer::vanLeer(const controller &cont, const runtimeMesh &mesh, flowModel &flow)
    : upwindScheme(cont, mesh, flow) {}
