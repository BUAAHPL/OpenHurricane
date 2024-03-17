/*!
 * \file HLLC.cpp
 * \brief Main subroutines for HLLC scheme.
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

#include "HLLC.hpp"
#include <cmath>

namespace OpenHurricane {
    createClassNameStr(HLLC, "HLLC");
    registerObjFty(upwindScheme, HLLC, controller);
} // namespace OpenHurricane

void OpenHurricane::HLLC::calcFlux(const real rhol, const real rhor, const vector &VL,
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

    const real &ul = VL[0];
    const real &ur = VR[0];
    const real &vl = VL[1];
    const real &vr = VR[1];
    const real &wl = VL[2];
    const real &wr = VR[2];

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

    // Total energy

    /*!\brief Total energy of left state.*/
    //real el = hl - pl / rhol;

    /*!\brief Total energy of right state.*/
    //real er = hr - pr / rhor;

    // Total enthalpy

    /*!\brief Total enthalpy of left state.*/
    real hl = el + pl / rhol;
    //real hl = hal + ekl;

    /*!\brief Total enthalpy of right state.*/
    real hr = er + pr / rhor;
    //real hr = har + ekr;

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
        qL = sqrt(1.0 + 0.5 * (gl + 1.0) * (ps / pl - 1.0) / gl);
    }
    if (ps > pr) {
        qR = sqrt(1.0 + 0.5 * (gr + 1.0) * (ps / pr - 1.0) / gr);
    }

    // Left wave speed SL
    real SL = unl - cl * qL;

    // Right wave speed SR
    real SR = unr + cr * qR;

    // Middle wave speed SM
    real sunl = SL - unl;
    real sunr = SR - unr;
    real SM = (pr - pl + rhol * unl * sunl - rhor * unr * sunr) / (rhol * sunl - rhor * sunr);

    real f = 0;
    real theta = 0;
    real MaL = 0;
    real MaR = 0;
    if (lowMachCor_ || enhancedShockStab_) {
        /*!\brief Mach number at left state.*/
        MaL = unl / cl;

        /*!\brief Mach number at right state.*/
        MaR = unr / cr;

        f = pow3(mag(min(pl / pr, pr / pl)));

        theta = min(max(mag(MaL), mag(MaR)), real(1));
    }

    /*if (flux.size() < 5)
    {
            flux.resize(5, Zero);
    }*/

    if (SL >= 0.0) {
        real rhou = rhol * unl;
        flux[0] = areaMag * rhou;
        flux[1] = areaMag * (rhou * ul + pl * nx);
        flux[2] = areaMag * (rhou * vl + pl * ny);
        flux[3] = areaMag * (rhou * wl + pl * nz);
        flux[4] = areaMag * rhou * hl;
    } else if (SL < 0.0 && SM >= 0.0) {
        //realArray ds(5, Zero);
        ////ds calculation
        //ds[1] = nx;
        //ds[2] = ny;
        //ds[3] = nz;
        //ds[4] = SM;

        //ds calculation
        ds_[0] = 0.0;
        ds_[1] = nx;
        ds_[2] = ny;
        ds_[3] = nz;
        ds_[4] = SM;

        real rhou = rhol * unl;
        flux[0] = rhou;
        flux[1] = (rhou * ul + pl * nx);
        flux[2] = (rhou * vl + pl * ny);
        flux[3] = (rhou * wl + pl * nz);
        flux[4] = rhou * hl;

        /*realArray uu(5, Zero);
        uu[0] = rhol;
        uu[1] = rhol * ul;
        uu[2] = rhol * vl;
        uu[3] = rhol * wl;
        uu[4] = rhol * el;*/

        uu_[0] = rhol;
        uu_[1] = rhol * ul;
        uu_[2] = rhol * vl;
        uu_[3] = rhol * wl;
        uu_[4] = rhol * el;

        real plr = pl + rhol * sunl * (SM - unl);
        if (lowMachCor_) {
            real plrss = theta * plr + (1.0 - theta) * (pl + pr) / 2.0;
            plr = f * plrss + (1.0 - f) * plr;
        }

        real mid = 1.0 / (SL - SM);
        for (integer i = 0; i < 5; ++i) {
            //real us = mid * (SL * uu[i] - flux[i] + plr * ds[i]);
            real us = mid * (SL * uu_[i] - flux[i] + plr * ds_[i]);
            flux[i] = areaMag * (flux[i] + SL * (us - uu_[i]));
        }
    } else if (SM < 0.0 && SR >= 0.0) {
        //realArray ds(5, Zero);
        ////ds calculation
        //ds[1] = nx;
        //ds[2] = ny;
        //ds[3] = nz;
        //ds[4] = SM;

        //ds calculation
        ds_[0] = 0.0;
        ds_[1] = nx;
        ds_[2] = ny;
        ds_[3] = nz;
        ds_[4] = SM;

        real rhou = rhor * unr;
        flux[0] = rhou;
        flux[1] = (rhou * ur + pr * nx);
        flux[2] = (rhou * vr + pr * ny);
        flux[3] = (rhou * wr + pr * nz);
        flux[4] = rhou * hr;

        /*realArray uu(5, Zero);
        uu[0] = rhor;
        uu[1] = rhor * ur;
        uu[2] = rhor * vr;
        uu[3] = rhor * wr;
        uu[4] = rhor * er;*/

        uu_[0] = rhor;
        uu_[1] = rhor * ur;
        uu_[2] = rhor * vr;
        uu_[3] = rhor * wr;
        uu_[4] = rhor * er;

        real plr = pr + rhor * sunr * (SM - unr);

        if (lowMachCor_) {
            real plrss = theta * plr + (1.0 - theta) * (pl + pr) / 2.0;
            plr = f * plrss + (1.0 - f) * plr;
        }

        real mid = 1.0 / (SR - SM);
        for (integer i = 0; i < 5; ++i) {
            //real us = mid * (SR * uu[i] - flux[i] + plr * ds[i]);
            real us = mid * (SR * uu_[i] - flux[i] + plr * ds_[i]);
            flux[i] = areaMag * (flux[i] + SR * (us - uu_[i]));
        }
    } else // SR < 0.0
    {
        real rhou = rhor * unr;
        flux[0] = areaMag * rhou;
        flux[1] = areaMag * (rhou * ur + pr * nx);
        flux[2] = areaMag * (rhou * vr + pr * ny);
        flux[3] = areaMag * (rhou * wr + pr * nz);
        flux[4] = areaMag * rhou * hr;
    }

    if (enhancedShockStab_) {
        real srhol = sqrt(rhol);
        real srhor = sqrt(rhor);

        real srholr = srhol + srhor;

        real cc = (srhol * cl + srhor * cr) / srholr;
        real uu = (srhol * ul + srhor * ur) / srholr;
        real vv = (srhol * vl + srhor * vr) / srholr;
        real ww = (srhol * wl + srhor * wr) / srholr;
        real MM = (srhol * MaL + srhor * MaR) / srholr;

        real phip = (f - 1.0) * SL * SR / (SR - SL) * inv(real(1) + mag(MM)) * (pr - pl) / sqr(cc);
        phip *= areaMag;

        flux[0] += phip;
        flux[1] += phip * uu;
        flux[2] += phip * vv;
        flux[3] += phip * ww;
        flux[4] += phip * 0.5 * (sqr(uu) + sqr(vv) + sqr(ww));
    }
}

OpenHurricane::HLLC::HLLC(const controller &cont, const runtimeMesh &mesh, flowModel &flow)
    : upwindScheme(cont, mesh, flow), ds_(5), uu_(5), lowMachCor_(false),
      enhancedShockStab_(false) {
    if (cont.found("HLLC")) {
        const auto &HLLCcont = cont.subController("HLLC");

        controllerSwitch myConts(HLLCcont);

        lowMachCor_ = myConts("lowMachCorrected", lowMachCor_);
        enhancedShockStab_ = myConts("enhancedShockStability", enhancedShockStab_);
    }
}
