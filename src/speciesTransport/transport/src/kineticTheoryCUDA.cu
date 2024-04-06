#include "kineticTheoryCUDA.hpp"
/*!
 * \file kineticTheoryCUDA.cu
 * \brief The subroutines and functions of kinetic theory transport in CUDA platform.
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

#include "kineticTheoryCUDA.hpp"
#include <cmath>
#ifdef CUDA_PARALLEL
#ifndef Rgen
#define Rgen cu_real(8.314e3)
#endif // !Rgen

#define muPre_ cu_real(2.6693e-6)
#define omegaPre_ cu_real(1.147)
#define omegaExponent_ cu_real(-0.145)

/**\brief Standard atmosphere pressure (Unit: [Pa]).*/
const cu_real OpenHurricane::CUDATransport::Patm = 1.01325e5;

/**\brief Universal gas constant (Unit: [J/(kmol K)]).*/
extern const cu_real OpenHurricane::CUDATransport::Ru = 8.3144626181532403e3;

cu_host OpenHurricane::CUDATransport::kineticTheoryCUDA::kineticTheoryCUDA(
    const cu_ushort nsp, const cu_real *__restrict__ ekb, const cu_real *__restrict__ sigma,
    const cuChem::speciesTable species)
    : ekb_(nsp, ekb), sigma_(nsp, sigma), species_(species) {}

cu_host OpenHurricane::CUDATransport::kineticTheoryCUDA::kineticTheoryCUDA(
    const cu_ushort nsp, const cu_real *__restrict__ ekb, const cu_real *__restrict__ sigma,
    const cuChem::speciesTable species, const cudaStreams &streams)
    : ekb_(nsp, ekb, streams()), sigma_(nsp, sigma, streams()), species_(species) {}

cu_device cu_real
OpenHurricane::CUDATransport::kineticTheoryCUDA::mu(const cu_real T, const cu_ushort isp) const {
    const cu_real Tstar = T / ekb_(isp);

    return muPre_ * sqrt(species_.Wi(isp) * T) / (cuSqr(sigma_(isp)) * omega22(Tstar, isp));
}

cu_device cu_real OpenHurricane::CUDATransport::kineticTheoryCUDA::Dim(
    const cu_real T, const cu_real p, const cu_real *__restrict__ xi,
    const cu_ushort isp) const {
    const auto ekbi = ekb_(isp);
    const auto sigmai = sigma_(isp);
    const auto Wi = species_.Wi(isp);
    const auto Tv1p5 = T * sqrt(T);

    cu_real XinvDij = 0;

    for (cu_ushort j = 0; j < nsp(); ++j) {
        const auto OMGDij = omegaDij(T, ekbi, ekb_(j));
        const cu_real sigmaij = 0.5 * (sigmai + sigma_(j));
        const cu_real Dij = cu_real(1.858e-7) * cu_real(1.01325e5) * Tv1p5 *
                              sqrt(cu_real(1) / Wi + cu_real(1) / species_.Wi(j)) /
                              (p * cuSqr(sigmaij) * OMGDij);

        const auto Xj = xi[j];
        if (isp != j) {
            XinvDij += (Xj / Dij);
        }
    }

    cu_real Dimm = 0;
    const auto Xii = xi[isp];
    if (Xii < cu_real(1) - cu_tiny) {
        Dimm = (cu_real(1) - Xii) / XinvDij;
    }
    if (isinf(Dimm)) {
        return 0;
    }
    return Dimm;
}

#endif // CUDA_PARALLEL