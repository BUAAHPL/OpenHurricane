/*!
 * \file kineticTheoryCUDA.hpp
 * \brief Header of kinetic theory transport in CUDA platform.
 *       The subroutines and functions are in the <i>kineticTheoryCUDA.cpp</i> file.
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

#include "speciesTableCUDA.hpp"

#include <cmath>
#ifndef Rgen
#define Rgen cu_real(8.3144626181532403e3)
#endif // !Rgen

namespace OpenHurricane {
    namespace CUDATransport {
        /**\brief Standard atmosphere pressure (Unit: [Pa]).*/
        extern const cu_real Patm;

        /**\brief Universal gas constant (Unit: [J/(kmol K)]).*/
        extern const cu_real Ru;

        /**
         * \brief CUDA kinetic theory transport.
         */
        class kineticTheoryCUDA {
        private:
            /*!
             *\brief An index indicating whether the molecule has a monatomic, linear or nonlinear geometrical configuration.
             *       If the index is 0, the molecule is a single atom. If the index is 1 the molecule is linear, and if it is 2, the molecule
             *       is nonlinear.
             */
            //cu1DArray<cu_integer> geometryIndex_;

            /*!\brief The Lennard-Jones potential well depth in Kelvins.*/
            cu1DArray<cu_real> ekb_;

            /*!\brief The Lennard-Jones collision diameter in angstroms.*/
            cu1DArray<cu_real> sigma_;

            /*!\brief The dipole moment in Debye.*/
            //cu1DArray<cu_real> mmu_;

            /*!\brief The polarizability in cubic angstroms.*/
            //cu1DArray<cu_real> alpha_;

            /*!\brief The rotational relaxation collision number at 298 K.*/
            //cu1DArray<cu_real> Zrot_;

            /** \brief The table of species. */
            const cuChem::speciesTable species_;

        public:
            /**
             * \brief Disallow null constructor.
             */
            inline cu_host kineticTheoryCUDA() = delete;

            cu_host kineticTheoryCUDA(const cu_ushort nsp,
                                             const cu_real *__restrict__ ekb,
                                             const cu_real *__restrict__ sigma,
                                             const cuChem::speciesTable species);

            cu_host kineticTheoryCUDA(const cu_ushort nsp,
                                             const cu_real *__restrict__ ekb,
                                             const cu_real *__restrict__ sigma,
                                             const cuChem::speciesTable species,
                                             const cudaStreams &streams);

            inline cu_dual kineticTheoryCUDA(const kineticTheoryCUDA &);

            /**
             * \brief Destructor.
             */
            cu_dual inline ~kineticTheoryCUDA() noexcept;

            /**
             * \brief Destroy.
             */
            inline cu_host void destroy();

            /*!
             * \brief Return the dynamic viscosity [kg/(m s)].
             * \param[in] T - Static temperature [K].
             * \param[in] isp - The index of species
             */
            cu_device cu_real mu(const cu_real T, const cu_ushort isp) const;

            /*!
             * \brief Return the thermal conductivity [W/mK].
             * \param[in] p - Static pressure [Pa].
             * \param[in] T - Static temperature [K].
             * \param[in] mui - Dynamic viscosity [kg/(m s)].
             * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
             */
            inline cu_device cu_real kappa(const cu_real T, const cu_real mui,
                                                    const cu_real cpi,
                                                    const cu_ushort isp) const;

            cu_device cu_real Dim(const cu_real T, const cu_real p,
                                           const cu_real *__restrict__ xi,
                                           const cu_ushort isp) const;

        private:
            /*!\brief The collision integral omega(2,2).*/
            inline cu_device cu_real omega22(const cu_real Tstar,
                                                      const cu_ushort isp) const;

            inline cu_device cu_real omegaDij(const cu_real T, const cu_real ekbi,
                                                       const cu_real ekbj) const;

        public:
            inline cu_device int nsp() const;

            /** \brief The table of species. */
            inline cu_device const cuChem::speciesTable &species() const noexcept;
        };
    } // namespace CUDATransport

    /**
     * \brief The class of reaction table.
     */
} // namespace OpenHurricane
#include "kineticTheoryCUDA.inl"

#ifdef Rgen
#undef Rgen
#endif // !Rgen

#endif // CUDA_PARALLEL
