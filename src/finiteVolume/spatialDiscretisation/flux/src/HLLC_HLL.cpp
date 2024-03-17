/*!
 * \file HLLC_HLL.cpp
 * \brief Main subroutines for HLLC_HLL scheme.
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

#include "HLLC_HLL.hpp"
#include "runtimeMesh.hpp"
#include <cmath>

namespace OpenHurricane {
    createClassNameStr(HLLC_HLL, "HLLC_HLL");
    registerObjFty(upwindScheme, HLLC_HLL, controller);
} // namespace OpenHurricane

void OpenHurricane::HLLC_HLL::calcFlux(const real rhol, const real rhor, const vector &VL,
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
    const real SL = unl - cl * qL;

    // Right wave speed SR
    const real SR = unr + cr * qR;

    // Middle wave speed SM
    const real sunl = SL - unl;
    const real sunr = SR - unr;
    const real SM = (pr - pl + rhol * unl * sunl - rhor * unr * sunr) / (rhol * sunl - rhor * sunr);

    //if (flux.size() < 5)
    //{
    //	flux.resize(5, Zero);
    //}

    if (SL >= 0.0) {
        const real rhou = rhol * unl;
        flux[0] = areaMag * rhou;
        flux[1] = areaMag * (rhou * ul + pl * nx);
        flux[2] = areaMag * (rhou * vl + pl * ny);
        flux[3] = areaMag * (rhou * wl + pl * nz);
        flux[4] = areaMag * rhou * hl;
    } else if (SL < 0.0 && SM >= 0.0) {
        //ds calculation
        ds_[0] = 0.0;
        ds_[1] = nx;
        ds_[2] = ny;
        ds_[3] = nz;
        ds_[4] = SM;

        real rhou = rhol * unl;
        fl_[0] = rhou;
        fl_[1] = (rhou * ul + pl * nx);
        fl_[2] = (rhou * vl + pl * ny);
        fl_[3] = (rhou * wl + pl * nz);
        fl_[4] = rhou * hl;

        uul_[0] = rhol;
        uul_[1] = rhol * ul;
        uul_[2] = rhol * vl;
        uul_[3] = rhol * wl;
        uul_[4] = rhol * el;

        rhou = rhor * unr;
        fr_[0] = rhou;
        fr_[1] = (rhou * ur + pr * nx);
        fr_[2] = (rhou * vr + pr * ny);
        fr_[3] = (rhou * wr + pr * nz);
        fr_[4] = rhou * hr;

        uur_[0] = rhor;
        uur_[1] = rhor * ur;
        uur_[2] = rhor * vr;
        uur_[3] = rhor * wr;
        uur_[4] = rhor * er;

        for (integer i = 0; i < 5; ++i) {
            uhll_[i] = (SR * uur_[i] - SL * uul_[i] + fl_[i] - fr_[i]) / (SR - SL);
        }
        const real plr = pl + rhol * sunl * (SM - unl);
        const real mid = 1.0 / (SL - SM);
        for (integer i = 0; i < 5; ++i) {
            uhllc_[i] = mid * (SL * uul_[i] - fl_[i] + plr * ds_[i]);
        }
        for (integer i = 0; i < 5; ++i) {
            uu_[i] = blend * uhllc_[i] + (1.0 - blend) * uhll_[i];
        }
        for (integer i = 0; i < 5; ++i) {
            flux[i] = areaMag * (fl_[i] + SL * (uu_[i] - uul_[i]));
        }
    } else if (SM < 0.0 && SR >= 0.0) {
        //ds calculation
        ds_[0] = 0.0;
        ds_[1] = nx;
        ds_[2] = ny;
        ds_[3] = nz;
        ds_[4] = SM;
        real rhou = rhol * unl;
        fl_[0] = rhou;
        fl_[1] = (rhou * ul + pl * nx);
        fl_[2] = (rhou * vl + pl * ny);
        fl_[3] = (rhou * wl + pl * nz);
        fl_[4] = rhou * hl;

        uul_[0] = rhol;
        uul_[1] = rhol * ul;
        uul_[2] = rhol * vl;
        uul_[3] = rhol * wl;
        uul_[4] = rhol * el;

        rhou = rhor * unr;
        fr_[0] = rhou;
        fr_[1] = (rhou * ur + pr * nx);
        fr_[2] = (rhou * vr + pr * ny);
        fr_[3] = (rhou * wr + pr * nz);
        fr_[4] = rhou * hr;

        uur_[0] = rhor;
        uur_[1] = rhor * ur;
        uur_[2] = rhor * vr;
        uur_[3] = rhor * wr;
        uur_[4] = rhor * er;
        for (integer i = 0; i < 5; ++i) {
            uhll_[i] = (SR * uur_[i] - SL * uul_[i] + fl_[i] - fr_[i]) / (SR - SL);
        }
        const real plr = pr + rhor * sunr * (SM - unr);
        const real mid = 1.0 / (SR - SM);
        for (integer i = 0; i < 5; ++i) {
            uhllc_[i] = mid * (SR * uur_[i] - fr_[i] + plr * ds_[i]);
        }
        for (integer i = 0; i < 5; ++i) {
            uu_[i] = blend * uhllc_[i] + (1.0 - blend) * uhll_[i];
        }
        for (integer i = 0; i < 5; ++i) {
            flux[i] = areaMag * (fr_[i] + SR * (uu_[i] - uur_[i]));
        }

    } else // SR < 0.0
    {
        const real rhou = rhor * unr;
        flux[0] = areaMag * rhou;
        flux[1] = areaMag * (rhou * ur + pr * nx);
        flux[2] = areaMag * (rhou * vr + pr * ny);
        flux[3] = areaMag * (rhou * wr + pr * nz);
        flux[4] = areaMag * rhou * hr;
    }
}
