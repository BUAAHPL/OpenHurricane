/*!
 * \file AUSMPlusUP.cpp
 * \brief Main subroutines for AUSM+ -up scheme.
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

#include "AUSMPlusUP.hpp"
#include <cmath>

namespace OpenHurricane {
    createClassNameStr(AUSMPlusUP, "AUSMPlusUP");
    registerObjFty(upwindScheme, AUSMPlusUP, controller);
} // namespace OpenHurricane
void OpenHurricane::AUSMPlusUP::calcFlux(const real rhol, const real rhor, const vector &VL,
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
    const real hl = el + pl / rhol;
    //real hl = hal + ekl;

    /*!\brief Total enthalpy of right state.*/
    const real hr = er + pr / rhor;
    //real hr = har + ekr;

    ///*!\brief Total enthalpy exclude chemical enthalpy of left state.*/
    //real htl = hl - hcl;

    ///*!\brief Total enthalpy exclude chemical enthalpy of right state.*/
    //real htr = hr - hcr;

    // Critical acoustic speed

    //// Critical acoustic speed of left state
    //real cs2 = 2.0 * (gl - 1.0) / (gl + 1.0) * htl;
    //real cs = std::sqrt(cs2);

    //real cla = cl * cl / OpenHurricane::max(cl, unl);

    //// Critical acoustic speed of left state
    //cs2 = 2.0 * (gr - 1.0) / (gr + 1.0) * htr;
    //cs = std::sqrt(cs2);

    //real cra = cr * cr / OpenHurricane::max(cr, -unr);

    //real cLR = OpenHurricane::min(cla, cra);
    const real cLR = 0.5 * (cl + cr);

    //// Total energy

    //	/*!\brief Total energy of left state.*/
    //	real el = hl - pl / rhol;

    //	/*!\brief Total energy of right state.*/
    //	real er = hr - pr / rhor;

    const real rhoLR = 0.5 * (rhol + rhor);

    // Mach number

    /*!\brief Mach number at left state.*/
    const real MaL = unl / cLR;

    /*!\brief Mach number at right state.*/
    const real MaR = unr / cLR;

    const real cLR2 = sqr(cLR);
    const real Maa2 = 0.5 * (sqr(unl) + sqr(unr)) / cLR2;
    const real Mo2 = min(real(1.0), max(Maa2, MaInf() * MaInf()));
    const real Mo = sqrt(Mo2);
    const real faMo = fa(Mo);

    const real Map =
        -Kp_ / faMo * max(real(1.0) - sigma_ * Maa2, Zero) * (pr - pl) / (rhoLR * cLR2);

    /*!\brief Mach number at the interface.*/
    const real ManLR = M4P(MaL) + M4M(MaR) + Map;

    // Pressure splitting

    const real _alpha = alpha(faMo);
    const real pu =
        -Ku_ * P5P(MaL, _alpha) * P5M(MaR, _alpha) * (rhol + rhor) * (faMo * cLR) * (unr - unl);

    /*!\brief Pressure flux at the interface.*/
    const real pLR = P5P(MaL, _alpha) * pl + P5M(MaR, _alpha) * pr + pu;

    /*!\brief Mass flux.*/
    real massLR;
    if (ManLR > 0.0) {
        massLR = cLR * ManLR * rhol;
    } else {
        massLR = cLR * ManLR * rhor;
    }

    /*!\brief Pressure flux.*/
    /*realArray fp(5, Zero);
    fp[1] = nx;
    fp[2] = ny;
    fp[3] = nz;
    fp = pLR * fp;*/
    fp_[0] = 0.0;
    fp_[1] = nx;
    fp_[2] = ny;
    fp_[3] = nz;
    fp_[4] = 0.0;
    fp_ = pLR * fp_;

    if (flux.size() < 5) {
        flux.resize(5, Zero);
    }
    /*!\brief Flux.*/
    if (massLR >= 0.0) {
        flux[0] = 1.0;
        flux[1] = ul;
        flux[2] = vl;
        flux[3] = wl;
        flux[4] = hl;

        flux *= massLR;
    } else {
        flux[0] = 1.0;
        flux[1] = ur;
        flux[2] = vr;
        flux[3] = wr;
        flux[4] = hr;

        flux *= massLR;
    }
    flux += fp_;
    flux *= areaMag;
}

OpenHurricane::AUSMPlusUP::AUSMPlusUP(const controller &cont, const runtimeMesh &mesh,
                                      flowModel &flow)
    : upwindScheme(cont, mesh, flow),
      Kp_(cont.subController("AUSMPlusUP").findOrDefault<real>("Kp", real(0.25))),
      Ku_(cont.subController("AUSMPlusUP").findOrDefault<real>("Ku", real(0.75))),
      sigma_(cont.subController("AUSMPlusUP").findOrDefault<real>("sigma", real(1.0))),
      beta_(cont.subController("AUSMPlusUP").findOrDefault<real>("beta", real(1.0 / 8.0))),
      fp_(5, Zero) {}
