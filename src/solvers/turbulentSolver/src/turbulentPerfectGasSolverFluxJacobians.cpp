/*!
 * \file turbulentPerfectGasSolverFluxJacobian.cpp
 * \brief Flux Jacobians of turbulent Perfect Gas Solver.
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

#include "turbulentPerfectGasSolver.hpp"

void OpenHurricane::turbulentPerfectGasSolver::Ac(const integer celli,
                                                  const vector &normal,
                                                  const real Vt,
                                                  realSquareMatrix &AC) const {
    solver::Ac(celli, normal, Vt, AC);
    if (turbPtr_->isCoupled()) {
        const real nx = normal.x();
        const real ny = normal.y();
        const real nz = normal.z();

        const real ux = v()[celli].x();
        const real vy = v()[celli].y();
        const real wz = v()[celli].z();

        const real V = normal * v()[celli];
        for (integer i = 0; i < turbPtr_->nEq(); ++i) {
            AC(i + 5, 0) = -V * turbPtr_->var(i)[celli];
            AC(i + 5, 1) = nx * turbPtr_->var(i)[celli];
            AC(i + 5, 2) = ny * turbPtr_->var(i)[celli];
            AC(i + 5, 3) = nz * turbPtr_->var(i)[celli];
            AC(i + 5, i + 5) = V;
        }
    }
}

void OpenHurricane::turbulentPerfectGasSolver::Acdq(const integer celli,
                                                    const vector &normal,
                                                    const realArray &dq,
                                                    realArray &adq) const {
    const real nx = normal.x();
    const real ny = normal.y();
    const real nz = normal.z();
    integer pqS = 5;
    if (turbPtr_->isCoupled()) {
        pqS += turbPtr_->nEq();
    }
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
    if (turbPtr_->isCoupled()) {
        for (integer i = 0; i < turbPtr_->nEq(); ++i) {
            adq[i + 5] = aa * turbPtr_->var(i)[celli] + uu * dq[i + 5];
        }
    }
}
