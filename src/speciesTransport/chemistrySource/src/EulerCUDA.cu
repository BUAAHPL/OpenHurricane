#include "EulerCUDA.hpp"
/*!
 * \file EulerCUDA.cpp
 * \brief Main subroutines for Euler ODE solver in CUDA.
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

#ifdef CUDA_PARALLEL
#include "EulerCUDA.hpp"

cu_device cu_real OpenHurricane::EulerCUDA::solve(
    CUDAChemistrySourceODEs &odes, const cu_real t0,
    const cu_real *__restrict__ y0, const cu_real *__restrict__ dydt0,
    const cu_real dt, cu_real *__restrict__ y,
    cu_real *__restrict__ y2tmp, cu_real *__restrict__ error,
    unsigned int tid, const unsigned int dimX) {
    auto j = tid;
    while (j < odes.nEqns()) {
        y2tmp[j] = y0[j] + 0.5 * dt * dydt0[j];
        //y[j] = y0[j] + dt * dydt0[j];
        j += dimX;
    }
    //return 0.5;
    __syncthreads();
    odes.DyDt(t0 + cu_real(0.5) * dt, y2tmp, error, y, tid, dimX);
    __syncthreads();

    j = tid;
    while (j < odes.nEqns()) {
        y2tmp[j] += 0.5 * dt * y[j];
        y[j] = y0[j] + dt * dydt0[j];
        cu_real yytmp = y2tmp[j];
        y2tmp[j] = y[j];
        y[j] = yytmp;
        error[j] = 2.0 / (1.0 - 2.0) * (y[j] - y2tmp[j]);
        j += dimX;
    }
    //__syncthreads();
    return maxError(odes, y2tmp, y, error, tid, dimX);
}

cu_device cu_real OpenHurricane::EulerCUDA::solveFactoring(
    CUDAChemistrySourceODEs &odes, const cu_real t0,
    const cu_real *__restrict__ y0, const cu_real *__restrict__ dydt0,
    const cu_real dt, const cu_real odeFactor, cu_real *__restrict__ y,
    cu_real *__restrict__ y2tmp, cu_real *__restrict__ error,
    unsigned int tid, const unsigned int dimX) {
    auto j = tid;
    while (j < odes.nEqns()) {
        //error[j] = dt * dydt0[j];
        //y[j] = y0[j] + error[j];
        y2tmp[j] = y0[j] + 0.5 * dt * odeFactor * dydt0[j];
        j += dimX;
    }
    __syncthreads();
    odes.DyDt(t0 + cu_real(0.5) * dt, y2tmp, error, y, tid, dimX);

    j = tid;
    while (j < odes.nEqns()) {
        y2tmp[j] += 0.5 * dt * odeFactor * y[j];
        y[j] = y0[j] + dt * odeFactor * dydt0[j];

        cu_real yytmp = y2tmp[j];
        y2tmp[j] = y[j];
        y[j] = yytmp;

        error[j] = 2.0 / (1.0 - 2.0) * (y[j] - y2tmp[j]);
        j += dimX;
    }
    __syncthreads();
    return maxError(odes, y2tmp, y, error, tid, dimX);
}

cu_device void OpenHurricane::EulerCUDA::solve(
    CUDAChemistrySourceODEs &odes, const cu_real t0,
    const cu_real *__restrict__ y0, const cu_real *__restrict__ dydt0,
    const cu_real dt, cu_real *__restrict__ y, unsigned int tid,
    const unsigned int dimX) {
    auto j = tid;
    while (j < odes.nEqns()) {
        y[j] = y0[j] + dt * dydt0[j];

        j += dimX;
    }
    __syncthreads();
}

cu_device void OpenHurricane::EulerCUDA::solveFactoring(
    CUDAChemistrySourceODEs &odes, const cu_real t0,
    const cu_real *__restrict__ y0, const cu_real *__restrict__ dydt0,
    const cu_real dt, const cu_real odeFactor, cu_real *__restrict__ y,
    unsigned int tid, const unsigned int dimX) {
    auto j = tid;
    while (j < odes.nEqns()) {
        y[j] = y0[j] + dt * odeFactor * dydt0[j];

        j += dimX;
    }
    __syncthreads();
}

cu_device void OpenHurricane::EulerCUDA::adaptiveSolve(
    CUDAChemistrySourceODEs &odes, cu_real &t, cu_real &dt0,
    const cu_real *__restrict__ y0, cu_real *__restrict__ y,
    cu_real *__restrict__ y2tmp, cu_real *__restrict__ dydt0,
    cu_real *__restrict__ error, unsigned int tid, const unsigned int dimX) {
    cu_real dt = dt0;
    odes.DyDt(t, y0, y, dydt0, tid, dimX);

    cu_real err = 0;
    do {
        err = solve(odes, t, y0, dydt0, dt, y, y2tmp, error, tid, dimX);

        __syncthreads();
        if (err > 1.0) {
            cu_real scale = cu_max(0.8 * pow(err, -0.25), 0.25);
            dt *= scale;
            if (dt <= cu_tiny) {
                dt = cu_tiny;
            }
        }
    } while (err > 1.0);
    t += dt;
    if (err > pow(cu_real(10) / cu_real(0.8),
                  cu_real(-1.0) / cu_real(0.25))) {
        dt0 = cu_min(
                  cu_max(cu_real(0.8) * pow(err, -cu_real(0.25)), 0.25),
                  10.0) *
              dt;
    } else {
        dt0 = cu_real(0.8) * cu_real(10) * dt;
    }
}

cu_device void OpenHurricane::EulerCUDA::adaptiveSolveFactoring(
    CUDAChemistrySourceODEs &odes, cu_real &t, cu_real &dt0,
    const cu_real odeFactor, const cu_real *__restrict__ y0,
    cu_real *__restrict__ y, cu_real *__restrict__ y2tmp,
    cu_real *__restrict__ dydt0, cu_real *__restrict__ error,
    unsigned int tid, const unsigned int dimX) {
    cu_real dt = dt0;
    odes.DyDt(t, y0, y, dydt0, tid, dimX);

    cu_real err = 0;
    do {
        err = solveFactoring(odes, t, y0, dydt0, dt, odeFactor, y, y2tmp, error,
                             tid, dimX);

        __syncthreads();
        if (err > 1.0) {
            cu_real scale = cu_max(0.8 * pow(err, -0.25), 0.25);
            dt *= scale;
            if (dt <= cu_tiny) {
                dt = cu_tiny;
            }
        }
    } while (err > 1.0);
    t += dt;
    if (err > pow(cu_real(10) / cu_real(0.8),
                  cu_real(-1.0) / cu_real(0.25))) {
        dt0 = cu_min(
                  cu_max(cu_real(0.8) * pow(err, -cu_real(0.25)), 0.25),
                  10.0) *
              dt;
    } else {
        dt0 = cu_real(0.8) * cu_real(10) * dt;
    }
}

cu_device void OpenHurricane::EulerCUDA::solve(
    CUDAChemistrySourceODEs &odes, const cu_real t0, cu_real &tn,
    cu_real &dt0, cu_real *__restrict__ y0, cu_real *__restrict__ y,
    cu_real *__restrict__ y2tmp, cu_real *__restrict__ dydt0,
    cu_real *__restrict__ error, unsigned int tid, const unsigned int dimX) {
    cu_real t = t0;
    cu_real dt = dt0;

    cu_real dtReached = t0;
    bool last = false;
    for (cu_integer n = 0; n < maxSteps_; ++n) {
        auto dtTry0 = dt;
        if ((t + dt - tn) * (t + dt - t0) > 0) {
            dt = tn - t;
            last = true;
        } else if (t + dt == tn) {
            last = true;
        }
        adaptiveSolve(odes, t, dt, y0, y, y2tmp, dydt0, error, tid, dimX);
        auto j = tid;
        while (j < odes.nEqns()) {
            y0[j] = y[j];
            j += dimX;
        }
        __syncthreads();
        dtReached = t;
        if ((t - tn) * (tn - t0) >= 0) {
            if (n > 0 && last) {
                dt0 = cu_min(dtTry0, tn);
            }
            dt = dt0;
            return;
        }
    }
    tn = dtReached;
    return;
}

cu_device void OpenHurricane::EulerCUDA::solveFactoring(
    CUDAChemistrySourceODEs &odes, const cu_real t0, cu_real &tn,
    cu_real &dt0, const cu_real odeFactor, cu_real *__restrict__ y0,
    cu_real *__restrict__ y, cu_real *__restrict__ y2tmp,
    cu_real *__restrict__ dydt0, cu_real *__restrict__ error,
    unsigned int tid, const unsigned int dimX) {
    cu_real t = t0;
    cu_real dt = dt0;

    cu_real dtReached = t0;
    bool last = false;
    for (cu_integer n = 0; n < maxSteps_; ++n) {
        auto dtTry0 = dt;
        if ((t + dt - tn) * (t + dt - t0) > 0) {
            dt = tn - t;
            last = true;
        }
        adaptiveSolveFactoring(odes, t, dt, odeFactor, y0, y, y2tmp, dydt0,
                               error, tid, dimX);
        auto j = tid;
        while (j < odes.nEqns()) {
            y0[j] = y[j];
            j += dimX;
        }
        __syncthreads();
        dtReached = t;
        if ((t - tn) * (tn - t0) >= 0) {
            if (n > 0 && last) {
                dt0 = dtTry0;
            }
            dt = dt0;
            return;
        }
    }
    tn = dtReached;
    return;
}

#endif // CUDA_PARALLEL