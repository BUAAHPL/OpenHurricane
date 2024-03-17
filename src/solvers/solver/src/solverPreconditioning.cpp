/*!
 * \file solverPreconditioning.cpp
 * \brief solverPreconditioning.
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

#include "solver.hpp"

void OpenHurricane::solver::PP(const integer cellI, realSquareMatrix &PM1,
                               realSquareMatrix &PP, realSquareMatrix &FM1,
                               realSquareMatrix &FP, real &aa4,
                               real &aa5) const {
    const auto &vel = v()[cellI];

    const real u = vel.x();
    const real v = vel.y();
    const real w = vel.z();

    // Derivative of the density with respect to the pressure
    const real rhop0 = mixtures().thermalTable().eos().DRhoDp(
        p()[cellI], rho()[cellI], T()[cellI], yi(), cellI);
    // Derivative of the density with respect to the temperature
    const real rhoT = mixtures().thermalTable().eos().DRhoDT(
        p()[cellI], rho()[cellI], T()[cellI], yi(), cellI);

    const real hp = (rho()[cellI] + rhoT * T()[cellI]) / sqr(rho()[cellI]);
    const real hT = mixtures().thermalTable().eos().cpm_p(
        p()[cellI], T()[cellI], yi(), cellI);
    const real H = E()[cellI] + p()[cellI] / rho()[cellI];
    const real q2 = sqr(u) + sqr(v) + sqr(w);

    auto pp = [&](const real rhop, realSquareMatrix &Pp,
                  realSquareMatrix &Pm) -> real {
        const real a1 =
            rho()[cellI] * rhop * hT + rhoT * (1.0 - rho()[cellI] * hp);

        // The transformation matrix from the conservative into the primitive variables
        Pp.setZero();
        {
            const real rhoTa1 = rhoT / a1;
            const real rhopa1 = rhop / a1;
            Pp(0, 0) = (rho()[cellI] * hT + rhoT * (H - q2)) / a1;
            Pp(0, 1) = rhoTa1 * u;
            Pp(0, 2) = rhoTa1 * v;
            Pp(0, 3) = rhoTa1 * w;
            Pp(0, 4) = -rhoTa1;
            Pp(1, 0) = -u / rho()[cellI];
            Pp(1, 1) = 1.0 / rho()[cellI];
            Pp(2, 0) = -v / rho()[cellI];
            Pp(2, 2) = 1.0 / rho()[cellI];
            Pp(3, 0) = -w / rho()[cellI];
            Pp(3, 3) = 1.0 / rho()[cellI];
            Pp(4, 0) = (1.0 - rhop * (H - q2) - rho()[cellI] * hp) / a1;
            Pp(4, 1) = -rhopa1 * u;
            Pp(4, 2) = -rhopa1 * v;
            Pp(4, 3) = -rhopa1 * w;
            Pp(4, 4) = rhopa1;
        }

        // The transformation matrix from the primitive into the conservative variables
        Pm.setZero();
        {
            Pm(0, 0) = rhop;
            Pm(0, 4) = rhoT;
            Pm(1, 0) = rhop * u;
            Pm(1, 1) = rho()[cellI];
            Pm(1, 4) = rhoT * u;
            Pm(2, 0) = rhop * v;
            Pm(2, 2) = rho()[cellI];
            Pm(2, 4) = rhoT * v;
            Pm(3, 0) = rhop * w;
            Pm(3, 3) = rho()[cellI];
            Pm(3, 4) = rhoT * w;
            Pm(4, 0) = rhop * H - 1 - rho()[cellI] * hp;
            Pm(4, 1) = rho()[cellI] * u;
            Pm(4, 2) = rho()[cellI] * v;
            Pm(4, 3) = rho()[cellI] * w;
            Pm(4, 4) = rhoT * H + rho()[cellI] * hT;
        }
        return a1;
    };

    const real aa1 = pp(rhop0, PM1, PP);

    const real delTah = pow(mesh().cellVolume()[cellI], real(1.0 / 3.0));
    const real eps = 1e-3;

    real dp = 0;

    for (integer fi = 0; fi < mesh().cells()[cellI].faceSize(); ++fi) {
        const auto &fac = mesh().faces()[mesh().cells()[cellI].facei(fi)];

        const integer lfc = fac.leftCell();
        const integer rgc = fac.rightCell();

        integer m = lfc;
        if (lfc == cellI) {
            m = rgc;
        } else {
            m = lfc;
        }

        dp += mag(p()[cellI] - p()[m]);
    }

    const real vv2 = sqrt(q2);
    real ur = max(max(vv2, max(mu()[cellI] / (rho()[cellI] * delTah),
                               kappal()[cellI] / delTah)),
                  eps * sqrt(dp / rho()[cellI]));
    const real ga = mixtures().thermalTable().eos().gammaThCorrectM(
        p()[cellI], rho()[cellI], T()[cellI], yi(), cellI);

    const real cc = ga * p()[cellI] / rho()[cellI];
    ur = min(ur, cc);

    const real theta = 1.0 / sqr(ur) - rhoT / (rho()[cellI] * hT);

    const real af1 = pp(theta, FM1, FP);

    aa4 = aa1 / af1;
    aa5 = rho()[cellI] * hT / af1;
}