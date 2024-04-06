/*!
 * \file JANAFCUDA.cu
 * \brief The subroutines and functions of JANAF thermo in CUDA platform.
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

#include "JANAFCUDA.hpp"
#include <cmath>
#ifdef CUDA_PARALLEL
#ifndef Rgen
#define Rgen cu_real(8.3144626181532403e3)
#endif // !Rgen
       /**\brief Standard atmosphere pressure (Unit: [Pa]).*/
const cu_real OpenHurricane::CUDAThermo::Patm = 1.01325e5;

/**\brief Universal gas constant (Unit: [J/(kmol K)]).*/
extern const cu_real OpenHurricane::CUDAThermo::Ru = 8.3144626181532403e3;

#endif // CUDA_PARALLEL

cu_host OpenHurricane::CUDAThermo::JANAFCUDA::JANAFCUDA(const cu_ushort nsp,
                                                           const cu_real *__restrict__ hhCoef,
                                                           const cu_real *__restrict__ hlCoef,
                                                           const cu_real *__restrict__ hTLow,
                                                           const cu_real *__restrict__ hTHigh,
                                                           const cu_real *__restrict__ hTComm)
    : highCpCoeffs_(nsp, 7, hhCoef), lowCpCoeffs_(nsp, 7, hlCoef), TLow_(nsp, hTLow),
      THigh_(nsp, hTHigh), TCommon_(nsp, hTComm) {}

OpenHurricane::CUDAThermo::JANAFCUDA::JANAFCUDA(const cu_ushort nsp,
                                            const cu_real *__restrict__ hhCoef,
                                            const cu_real *__restrict__ hlCoef,
                                            const cu_real *__restrict__ hTLow,
                                            const cu_real *__restrict__ hTHigh,
                                            const cu_real *__restrict__ hTComm,
                                            const cudaStreams &streams)
    : highCpCoeffs_(nsp, 7), lowCpCoeffs_(nsp, 7), TLow_(nsp, hTLow, streams()),
      THigh_(nsp, hTHigh, streams()), TCommon_(nsp, hTComm, streams()) {
    highCpCoeffs_.copyFromHostAsync(hhCoef, streams());
    lowCpCoeffs_.copyFromHostAsync(hlCoef, streams());
}

cu_device cu_real OpenHurricane::CUDAThermo::JANAFCUDA::Cp0(const cu_real T,
                                                                 const cu_ushort isp) const {
    const auto &a = coeff(T, isp);
    const auto T0 = limit(T, isp);
    return Rgen *
           ((((a(4, isp) * T0 + a(3, isp)) * T0 + a(2, isp)) * T0 + a(1, isp)) * T0 + a(0, isp));
}

cu_device cu_real OpenHurricane::CUDAThermo::JANAFCUDA::Ha0(const cu_real T,
                                                                 const cu_ushort isp) const {
    const auto &a = coeff(T, isp);
    if (T < TLow_(isp)) {
        const cu_real T0 = TLow_(isp);

        cu_real H0 =
            Rgen * (((((a(4, isp) / 5.0 * T0 + a(3, isp) / 4.0) * T0 + a(2, isp) / 3.0) * T0 +
                      a(1, isp) / 2.0) *
                         T0 +
                     a(0, isp)) *
                        T0 +
                    a(5, isp));
        const cu_real Cp0 =
            Rgen *
            ((((a(4, isp) * T0 + a(3, isp)) * T0 + a(2, isp)) * T0 + a(1, isp)) * T0 + a(0, isp));
        H0 += Cp0 * (T - T0);

        return H0;
    } else if (T > THigh_(isp)) {
        const cu_real T0 = THigh_(isp);

        cu_real H0 =
            Rgen * (((((a(4, isp) / 5.0 * T0 + a(3, isp) / 4.0) * T0 + a(2, isp) / 3.0) * T0 +
                      a(1, isp) / 2.0) *
                         T0 +
                     a(0, isp)) *
                        T0 +
                    a(5, isp));
        const cu_real Cp0 =
            Rgen *
            ((((a(4, isp) * T0 + a(3, isp)) * T0 + a(2, isp)) * T0 + a(1, isp)) * T0 + a(0, isp));
        H0 += Cp0 * (T - T0);

        return H0;
    } else {
        return Rgen * (((((a(4, isp) / 5.0 * T + a(3, isp) / 4.0) * T + a(2, isp) / 3.0) * T +
                         a(1, isp) / 2.0) *
                            T +
                        a(0, isp)) *
                           T +
                       a(5, isp));
    }
}
