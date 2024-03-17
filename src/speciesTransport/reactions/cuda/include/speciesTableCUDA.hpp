/*!
 * \file speciesTableCUDA.hpp
 * \brief Header of species list in CUDA platform.
 *       The subroutines and functions are in the <i>speciesTableCUDA.cpp</i> file.
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
#pragma once

#ifdef CUDA_PARALLEL
#include "CUDAFunctions.hpp"
#include "JANAFCUDA.hpp"
#include "physicalThermoConatsnts.hpp"
#include <cmath>

#include "cudaStreams.hpp"
#include "cuChemInterface.hpp"

#ifndef Rgen
#define Rgen cu_real(8.3144626181532403e3)
#endif // !Rgen
namespace OpenHurricane {
namespace cuChem {
class speciesTable {
private:
    /** \brief Molercular weight array with size of nsp. */
    cu1DArray<cu_real> W_;
    CUDAThermo::JANAFCUDA thermo_;

public:
    /**
     * \brief Null constructor.
     */
    inline cu_host speciesTable();

    /**
     * \brief Construct from components.
     * \param[in] nsp - Number of species
     * \param[in] WPtr - The pointer of molercular weight array with size of nrc
     * \param[in] TCommonPtr - The pointer of array with size of nrcPD
     * \param[in] TLowPtr - The pointer of array with size of nrcPD
     * \param[in] THighPtr - The pointer of array with size of nrcPD
     * \param[in] highCpCoeffsPtr - The pointer of 2D array with size of (7,nsp)
     * \param[in] lowCpCoeffsPtr - The pointer of 2D array with size of (7,nsp)
     */
    inline cu_host
    speciesTable(const cu_ushort nsp, const cu_real *__restrict__ WPtr,
                 const cu_real *__restrict__ TCommonPtr,
                 const cu_real *__restrict__ TLowPtr,
                 const cu_real *__restrict__ THighPtr,
                 const cu_real *__restrict__ highCpCoeffsPtr,
                 const cu_real *__restrict__ lowCpCoeffsPtr);

    /**
     * \brief Construct from components.
     * \param[in] nsp - Number of species
     * \param[in] WPtr - The pointer of molercular weight array with size of nrc
     * \param[in] TCommonPtr - The pointer of array with size of nrcPD
     * \param[in] TLowPtr - The pointer of array with size of nrcPD
     * \param[in] THighPtr - The pointer of array with size of nrcPD
     * \param[in] highCpCoeffsPtr - The pointer of 2D array with size of (7,nsp)
     * \param[in] lowCpCoeffsPtr - The pointer of 2D array with size of (7,nsp)
     */
    inline cu_host
    speciesTable(const cu_ushort nsp, const cu_real *__restrict__ WPtr,
                 const cu_real *__restrict__ TCommonPtr,
                 const cu_real *__restrict__ TLowPtr,
                 const cu_real *__restrict__ THighPtr,
                 const cu_real *__restrict__ highCpCoeffsPtr,
                 const cu_real *__restrict__ lowCpCoeffsPtr,
                 const cudaStreams &streams);

    /**
     * \brief Construct from reactionCUDAInterface::speciesTableInterface.
     */
    inline cu_host
    speciesTable(const cu_ushort nsp,
                 const cuChemInterface::speciesTableInterface &sptInt);

    /**
     * \brief Construct from reactionCUDAInterface::speciesTableInterface.
     */
    inline cu_host
    speciesTable(const cu_ushort nsp,
                 const cuChemInterface::speciesTableInterface &sptInt,
                 const cudaStreams &streams);

    /**
     * \brief Copy constructor.
     */
    inline cu_dual speciesTable(const speciesTable &spt);

    /**
     * \brief Destructor.
     */
    cu_dual inline ~speciesTable() noexcept;

    /**
     * \brief Destroy cuda array memory.
     */
    inline cu_host void destroy();

    /**
     * \brief Calculate negative Gibbâ€™s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
     */
    inline cu_device cu_real nGdRT(const cu_ushort isp,
                                            const cu_real Td) const;

    inline cu_dual void yiToci(const cu_real rho,
                                      const cu_real *__restrict__ yi,
                                      const cu_real *__restrict__ wi,
                                      const cu_ushort isp,
                                      cu_real *__restrict__ ci) const;

    inline cu_device void yiToci(const cu_real rho,
                                        const cu_real *__restrict__ yi,
                                        const cu_ushort isp,
                                        cu_real *__restrict__ ci) const;

    /*inline cu_dual cu_real yiToci
    (
            const cu_real rho,
            const cu_real yi,
            const cu_real wi
    )const;*/

    inline cu_device cu_real yiToci(const cu_real rho,
                                             const cu_real yi,
                                             const cu_ushort isp) const;

    inline cu_device cu_real Wi(const cu_ushort isp) const;

    inline cu_device const CUDAThermo::JANAFCUDA &thermo() const;

    /*!\brief Heat capacity at constant pressure [J/(kg K)] at standard pressure (1 atm).*/
    inline cu_device cu_real cp0(const cu_real T,
                                          const cu_ushort isp) const;

    /*!\brief Molar heat capacity at constant pressure [J/(kmol K)] at standard pressure (1 atm).*/
    inline cu_device cu_real Cp0(const cu_real T,
                                          const cu_ushort isp) const;

    /*!\brief Heat capacity at constant volume [J/(kg K)] at standard pressure (1 atm).*/
    inline cu_device cu_real cv0(const cu_real T,
                                          const cu_ushort isp) const;

    /*!\brief Molar heat capacity at constant volume [J/(kmol K)] at standard pressure (1 atm).*/
    inline cu_device cu_real Cv0(const cu_real T,
                                          const cu_ushort isp) const;

    /*!\brief Absolute Enthalpy [J/kg] at standard pressure.*/
    inline cu_device cu_real ha0(const cu_real T,
                                          const cu_ushort isp) const;

    /*!\brief Molar absolute Enthalpy [J/kmol] at standard pressure.*/
    inline cu_device cu_real Ha0(const cu_real T,
                                          const cu_ushort isp) const;

    /*!\brief Absolute internal energy [J/kg] at standard pressure.*/
    inline cu_device cu_real ea0(const cu_real T,
                                          const cu_ushort isp) const;

    /*!\brief Molar absolute internal energy [J/kmol] at standard pressure.*/
    inline cu_device cu_real Ea0(const cu_real T,
                                          const cu_ushort isp) const;

    /*!\brief Entropy [J/(kg K)] at standard pressure.*/
    inline cu_device cu_real s0(const cu_real T,
                                         const cu_ushort isp) const;

    /*!\brief Molar entropy [J/(kmol K)] at standard pressure.*/
    inline cu_device cu_real S0(const cu_real T,
                                         const cu_ushort isp) const;

    /** \brief Gas constant [J/(kg K)]*/
    inline cu_device cu_real Ri(const cu_ushort isp) const;

    inline cu_device void yidWi(const cu_ushort isp,
                                       const cu_real *__restrict__ yi,
                                       cu_real *__restrict__ ydwi) const;

    inline cu_device void yi2xi(const cu_ushort isp,
                                       const cu_real *__restrict__ yi,
                                       const cu_real wm,
                                       cu_real *__restrict__ xi) const;
};

} // namespace CUDAReactions
} // namespace OpenHurricane

#include "speciesTableCUDA.inl"

#ifdef Rgen
#undef Rgen
#endif // !Rgen
#endif // CUDA_PARALLEL
