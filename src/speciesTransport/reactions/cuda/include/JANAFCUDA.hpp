/*!
 * \file JANAFCUDA.hpp
 * \brief Header of JANAF thermo in CUDA platform.
 *       The subroutines and functions are in the <i>JANAFCUDA.cpp</i> file.
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
#pragma once

#ifdef CUDA_PARALLEL
#include "CUDAFunctions.hpp"
#include "cudaStreams.hpp"
#include <cmath>
#ifndef Rgen
#define Rgen cu_real(8.3144626181532403e3)
#endif // !Rgen

namespace OpenHurricane {
    namespace CUDAThermo {
        /**\brief Standard atmosphere pressure (Unit: [Pa]).*/
        extern const cu_real Patm;

        /**\brief Universal gas constant (Unit: [J/(kmol K)]).*/
        extern const cu_real Ru;

        /**
         * \brief CUDA JANAF thermo.
         */
        class JANAFCUDA {
        private:
            cu2DArray<cu_real> highCpCoeffs_;
            cu2DArray<cu_real> lowCpCoeffs_;

            cu1DArray<cu_real> TLow_;
            cu1DArray<cu_real> THigh_;
            cu1DArray<cu_real> TCommon_;

        public:
            /**
             * \brief Disallow null constructor.
             */
            inline cu_host JANAFCUDA();

            cu_host JANAFCUDA(const cu_ushort nsp, const cu_real *__restrict__ hhCoef,
                                     const cu_real *__restrict__ hlCoef,
                                     const cu_real *__restrict__ hTLow,
                                     const cu_real *__restrict__ hTHigh,
                                     const cu_real *__restrict__ hTComm);

            cu_host JANAFCUDA(const cu_ushort nsp, const cu_real *__restrict__ hhCoef,
                                     const cu_real *__restrict__ hlCoef,
                                     const cu_real *__restrict__ hTLow,
                                     const cu_real *__restrict__ hTHigh,
                                     const cu_real *__restrict__ hTComm,
                                     const cudaStreams &streams);

            inline cu_dual JANAFCUDA(const JANAFCUDA &);

            cu_dual inline ~JANAFCUDA() noexcept;

            inline cu_host void destroy();

            /**
             * \brief Molar heat capacity at constant pressure [J/(kmol K)] at standard pressure (1 atm).
             * \param[in] T - Temperature [K]
             * \param[in] isp - Index of species
             */
            cu_device cu_real Cp0(const cu_real T, const cu_ushort isp) const;

            /**
             * \brief Molar absolute Enthalpy [J/kmol] at standard pressure.
             * \param[in] T - Temperature
             * \param[in] isp - Index of species
             */
            cu_device cu_real Ha0(const cu_real T, const cu_ushort isp) const;

            /**
             * \brief Molar entropy [J/(kmol K)] at standard pressure.
             */
            inline cu_device cu_real S0(const cu_real T, const cu_ushort isp) const;

            /**
             * \brief Return \f${\frac{{TS_i^0 - H_i^0}}{{{R_0}T}}}\f$.
             *        i.e. \f$-{\frac{{G_i^0}}{{{R_0}T}}}\]\f$
             */
            cu_device cu_real SH(const cu_real T, const cu_ushort isp) const;

            /**
             * \brief Origin coefficients from JANAF type thermo file.
             */
            inline cu_device const cu2DArray<cu_real> &
            coeff(const cu_real T, const cu_ushort isp) const;

            inline cu_device cu_real limit(const cu_real T, const cu_ushort isp) const;
        };
    } // namespace CUDAThermo

    /**
     * \brief The class of reaction table.
     */
} // namespace OpenHurricane
#include "JANAFCUDA.inl"

#ifdef Rgen
#undef Rgen
#endif // !Rgen

#endif // CUDA_PARALLEL
