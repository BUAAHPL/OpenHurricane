/*!
 * \file transportListCUDA.hpp
 * \brief Header of transport in CUDA platform.
 *       The subroutines and functions are in the <i>transportListCUDA.cpp</i> file.
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
#include "cudaStreams.hpp"
#include "kineticTheoryCUDA.hpp"
#include "nThreadsAndBlocks.hpp"

namespace OpenHurricane {
    namespace CUDATransport {
        /**
         * \brief CUDA JANAF thermo.
         */
        class transportListCUDA {
        private:
            kineticTheoryCUDA transport_;

        public:
            /**
             * \brief Null constructor.
             */
            inline cu_host transportListCUDA() = delete;

            cu_host transportListCUDA(const cu_ushort nsp,
                                              const cu_real *__restrict__ ekb,
                                              const cu_real *__restrict__ sigma,
                                              const cuChem::speciesTable species);

            cu_host transportListCUDA(const cu_ushort nsp,
                                              const cu_real *__restrict__ ekb,
                                              const cu_real *__restrict__ sigma,
                                              const cuChem::speciesTable species,
                                              const cudaStreams &streams);

            inline cu_dual transportListCUDA(const transportListCUDA &);

            /**
             * \brief Destructor.
             */
            cu_dual inline ~transportListCUDA() noexcept;

            /**
             * \brief Destroy.
             */
            inline cu_host void destroy();

            inline cu_dual const kineticTheoryCUDA &transport() const noexcept;

            /**
             * \brief Get molercular viscosity of mixtures by Wilke's semi-empirical formula.
             */
            inline cu_device void mu(cu_real *__restrict__ mui,
                                            const cu_real *__restrict__ xi,
                                            const cu_real *__restrict__ wci,
                                            const cu_ushort i) const;

            /**
             * \brief Get molercular viscosity of mixtures by Wilke's semi-empirical formula.
             */
            inline cu_device void muKappa(cu_real *__restrict__ mui,
                                                 cu_real *__restrict__ kappai,
                                                 const cu_real *__restrict__ xi,
                                                 const cu_real *__restrict__ wci,
                                                 const cu_ushort i) const;

            cu_device void WilkeCoeff(const cu_real *__restrict__ mui,
                                             const cu_real *__restrict__ xi, const cu_ushort i,
                                             cu_real *__restrict__ wci) const;

        private:
            inline cu_device cu_real PhiMuij(const cu_real mui, const cu_real muj,
                                                      const cu_real Wi, const cu_real Wj) const;

            inline cu_device cu_real PhiKij(const cu_real mui, const cu_real muj,
                                                     const cu_real Wi, const cu_real Wj) const;

        public:
            inline cu_device int nsp() const;

            /** \brief The table of species. */
            inline cu_device const cuChem::speciesTable &species() const noexcept;
        };

        void calcTransportProperties(const cu_real *__restrict__ hostYiTP,
                                     const transportListCUDA &tran, const nThreadsAndBlocks &nTB,
                                     const cu_integer nsp, const cu_integer nCells,
                                     cu_real *__restrict__ hostDimMuKappa);

        /**
         * \brief Get transport properties.
         * \note yiT would be changed into xiT after being called.
         */
        void calcTransportPropertiesAsync(const cu_real *__restrict__ hostYiTP,
                                          cu2DArray<cu_real> &dYiTP,
                                          const transportListCUDA &tran,
                                          const nThreadsAndBlocks &nTB, const cu_integer nsp,
                                          const cu_integer nCells,
                                          cu_real *__restrict__ hostDimMuKappa,
                                          cu2DArray<cu_real> &dDimMuKappa,
                                          const cudaStreams &streams, const bool copyDYiT = false);
    } // namespace CUDATransport

} // namespace OpenHurricane
#include "transportListCUDA.inl"

#ifdef Rgen
#undef Rgen
#endif // !Rgen

#endif // CUDA_PARALLEL
