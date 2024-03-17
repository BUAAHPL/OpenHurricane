/*!
 * \file spatialSchemeLowMachCorr.cpp
 * \brief Main subroutines for spatial scheme with Dimitris's Low-Mach number correction: change velocity, keep
          energy constant & change static pressure in subsonic flows.
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

#include "spatialScheme.hpp"

void OpenHurricane::spatialScheme::lowMachCorrection(const real rhol, const real rhor,
                                                     const real gammal, const real gammar,
                                                     vector &vl, vector &vr, real &pl, real &pr,
                                                     const vector &fA) const {
    const real gm1l = gammal - 1;
    const real gm1r = gammar - 1;

    //total energy, not inculd chemical energy
    const real eel = pl / gm1l + real(0.5) * rhol * vl.magSqr();
    const real eer = pr / gm1r + real(0.5) * rhor * vr.magSqr();

    //face normal Mach number
    const auto n = fA.normalized();
    const real Manl = mag(vl * n) / sqrt(gammal * pl / rhol);
    const real Manr = mag(vr * n) / sqrt(gammar * pr / rhor);
    const real Manf = max(Manl, Manr);
    if (Manf >= real(1) || Manf == real(0)) {
        return;
    }

    //velocity correction for subsonic flows
    const real zeta = min(real(1), Manf);
    const real OpZ = real(1.0) + zeta;
    const real OmZ = real(1.0) - zeta;
    const auto vlc = real(0.50) * (OpZ * vl + OmZ * vr);
    const auto vrc = real(0.50) * (OpZ * vr + OmZ * vl);
    vl = vlc;
    vr = vrc;

    //keep energy constant and change static pressure
    pl = gm1l * (eel - 0.5 * rhol * vl.magSqr());
    pr = gm1r * (eer - 0.5 * rhor * vr.magSqr());
}
