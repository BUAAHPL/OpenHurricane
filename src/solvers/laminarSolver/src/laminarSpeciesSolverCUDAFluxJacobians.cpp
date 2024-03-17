/*!
 * \file laminarSpeciesSolverCUDA.cpp
 * \brief The flux Jacobians of laminar species transport Solver using CUDA accelerating chemical source terms calculation.
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

#include "laminarSpeciesSolverCUDA.hpp"
#ifdef CUDA_PARALLEL
void OpenHurricane::laminarSpeciesSolverCUDA::Ac(const integer celli,
                                                 const vector &normal,
                                                 const real Vt,
                                                 realSquareMatrix &AC) const {
    AC.setZero();
    const integer nsp = specTable().size();

    integer pqS = 5;
    pqS += nsp - 1;
    realArray pq(pqS);
    dpdq(celli, normal, pq);

    const real nx = normal.x();
    const real ny = normal.y();
    const real nz = normal.z();

    const real ux = v()[celli].x();
    const real vy = v()[celli].y();
    const real wz = v()[celli].z();

    const real h = E()[celli] + p()[celli] / rho()[celli];
    const real V = normal * v()[celli];

    {
        AC(0, 0) = 0;
        AC(0, 1) = nx;
        AC(0, 2) = ny;
        AC(0, 3) = nz;
        AC(0, 4) = 0.0;
        AC(1, 0) = pq[0] * nx - ux * V;
        AC(1, 1) = pq[1] * nx + V + ux * nx;
        AC(1, 2) = pq[2] * nx + ux * ny;
        AC(1, 3) = pq[3] * nx + ux * nz;
        AC(1, 4) = pq[4] * nx;
        AC(2, 0) = pq[0] * ny - vy * V;
        AC(2, 1) = pq[1] * ny + vy * nx;
        AC(2, 2) = pq[2] * ny + V + vy * ny;
        AC(2, 3) = pq[3] * ny + vy * nz;
        AC(2, 4) = pq[4] * ny;
        AC(3, 0) = pq[0] * nz - wz * V;
        AC(3, 1) = pq[1] * nz + wz * nx;
        AC(3, 2) = pq[2] * nz + wz * ny;
        AC(3, 3) = pq[3] * nz + V + wz * nz;
        AC(3, 4) = pq[4] * nz;
        AC(4, 0) = pq[0] * V - h * V;
        AC(4, 1) = pq[1] * V + h * nx;
        AC(4, 2) = pq[2] * V + h * ny;
        AC(4, 3) = pq[3] * V + h * nz;
        AC(4, 4) = V + pq[4] * V;
    }

    integer offset = 5;

    for (integer i = 0; i < nsp - 1; ++i) {
        AC(i + offset, 0) = -V * yi()[i][celli];
        AC(i + offset, 1) = nx * yi()[i][celli];
        AC(i + offset, 2) = ny * yi()[i][celli];
        AC(i + offset, 3) = nz * yi()[i][celli];

        AC(1, i + offset) = nx * pq[i + offset];
        AC(2, i + offset) = ny * pq[i + offset];
        AC(3, i + offset) = nz * pq[i + offset];
        AC(4, i + offset) = V * pq[i + offset];

        AC(i + offset, i + offset) = V;
    }
}

void OpenHurricane::laminarSpeciesSolverCUDA::dpdq(const integer celli,
                                                   const vector &normal,
                                                   realArray &pq) const {
    solver::dpdq(celli, normal, pq);
    const real gm1 = gama()[celli] - 1.0;
    const integer nsp = specTable().size();
    const real RinspT = specTable().Ri(nsp - 1) * T()[celli];
    const real ensp = hi()[nsp - 1][celli] - RinspT;
    pq[0] += (RinspT - gm1 * ensp);
    for (integer i = 0; i < nsp - 1; ++i) {
        const real RiiT = specTable().Ri(i) * T()[celli];
        const real ei = hi()[i][celli] - RiiT;
        pq[5 + i] = (RiiT - RinspT) - gm1 * (ei - ensp);
    }
}

void OpenHurricane::laminarSpeciesSolverCUDA::Acdq(const integer celli,
                                                   const vector &normal,
                                                   const realArray &dq,
                                                   realArray &adq) const {
    const real nx = normal.x();
    const real ny = normal.y();
    const real nz = normal.z();
    const integer nsp = specTable().size();
    integer pqS = 5;
    pqS += nsp - 1;
    realArray pq(pqS);
    dpdq(celli, normal, pq);
    const real ux = v()[celli].x();
    const real vy = v()[celli].y();
    const real wz = v()[celli].z();
    const real ek = 0.5 * v()[celli].magSqr();
    const real uu = ux * nx + vy * ny + wz * nz;
    const real aa = nx * dq[1] + ny * dq[2] + nz * dq[3] - uu * dq[0];
    real bb = Zero;
    for (integer i = 0; i < pq.size(); ++i) {
        bb += dq[i] * pq[i];
    }

    const real h = E()[celli] + p()[celli] / rho()[celli];

    adq[0] = aa + uu * dq[0];
    adq[1] = aa * ux + uu * dq[1] + bb * nx;
    adq[2] = aa * vy + uu * dq[2] + bb * ny;
    adq[3] = aa * wz + uu * dq[3] + bb * nz;
    adq[4] = aa * h + uu * dq[4] + bb * uu;

    for (integer i = 0; i < nsp - 1; ++i) {
        adq[i + 5] = aa * yi()[i][celli] + uu * dq[i + 5];
    }
}

#endif // CUDA_PARALLEL