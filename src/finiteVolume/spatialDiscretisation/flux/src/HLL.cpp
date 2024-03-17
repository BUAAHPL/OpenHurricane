/*!
 * \file HLL.cpp
 * \brief Main subroutines for HLL scheme.
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

#include "HLL.hpp"
#include <cmath>

namespace OpenHurricane {
    createClassNameStr(HLL, "HLL");
    registerObjFty(upwindScheme, HLL, controller);
} // namespace OpenHurricane

void OpenHurricane::HLL::calcFlux(const real rhol, const real rhor, const vector &VL,
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

    // Total energy

    /*!\brief Total energy of left state.*/
    //real el = hl - pl / rhol;

    /*!\brief Total energy of right state.*/
    //real er = hr - pr / rhor;

    // The density of the interface
    real rhoAve = 0.5 * (rhol + rhor);

    // The acoustic speed of the interface
    real cAve = 0.5 * (cl + cr);

    // Estimate p*
    real ps = 0.5 * (pl + pr) - 0.5 * (unr - unl) * rhoAve * cAve;
    ps = max(Zero, ps);

    // Estimate wave speeds from p*

    real qL = 1.0;
    real qR = 1.0;
    if (ps > pl) {
        qL = std::sqrt(1.0 + 0.5 * (gl + 1.0) * (ps / pl - 1.0) / gl);
    }
    if (ps > pr) {
        qR = std::sqrt(1.0 + 0.5 * (gr + 1.0) * (ps / pr - 1.0) / gr);
    }

    // Left wave speed SL
    real SL = unl - cl * qL;

    // Right wave speed SR
    real SR = unr + cr * qR;

    if (flux.size() < 5) {
        flux.resize(5, Zero);
    }

    if (SL >= 0.0) {
        real rhou = rhol * unl;
        flux[0] = areaMag * rhou;
        flux[1] = areaMag * (rhou * ul + pl * nx);
        flux[2] = areaMag * (rhou * vl + pl * ny);
        flux[3] = areaMag * (rhou * wl + pl * nz);
        flux[4] = areaMag * rhou * hl;
    } else if (SL < 0.0 && SR > 0.0) {
        real rhou = rhol * unl;
        flux[0] = rhou;
        flux[1] = (rhou * ul + pl * nx);
        flux[2] = (rhou * vl + pl * ny);
        flux[3] = (rhou * wl + pl * nz);
        flux[4] = rhou * hl;

        rhou = rhor * unr;

        flux[0] = SR * flux[0] - SL * rhou;
        flux[1] = SR * flux[1] - SL * (rhou * ur + pr * nx);
        flux[2] = SR * flux[2] - SL * (rhou * vr + pr * ny);
        flux[3] = SR * flux[3] - SL * (rhou * wr + pr * nz);
        flux[4] = SR * flux[4] - SL * rhou * hr;

        realArray uu(5, Zero);

        uu[0] = SR * SL * (rhor - rhol);
        uu[1] = SR * SL * (rhor * ur - rhol * ul);
        uu[2] = SR * SL * (rhor * vr - rhol * vl);
        uu[3] = SR * SL * (rhor * wr - rhol * wl);
        uu[4] = SR * SL * (rhor * er - rhol * el);

        flux = (areaMag / (SR - SL)) * (flux + uu);

    } else // SR <= 0.0
    {
        real rhou = rhor * unr;
        flux[0] = areaMag * rhou;
        flux[1] = areaMag * (rhou * ur + pr * nx);
        flux[2] = areaMag * (rhou * vr + pr * ny);
        flux[3] = areaMag * (rhou * wr + pr * nz);
        flux[4] = areaMag * rhou * hr;
    }
}
