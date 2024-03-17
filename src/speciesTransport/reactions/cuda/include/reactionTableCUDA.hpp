/*!
 * \file reactionTableCUDA.hpp
 * \brief Header of reaction table in CUDA platform.
 *       The subroutines and functions are in the <i>reactionTableCUDA.cpp</i> file.
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
#include "cuChemInterface.hpp"
#include "cudaStreams.hpp"
#include "physicalThermoConatsnts.hpp"
#include "reactionRateCoeffsCUDA.hpp"
#include "reactionRatePressureDependent.hpp"
#include "reactionRateThirdBody.hpp"
#include "speciesTableCUDA.hpp"
#include <cmath>

namespace OpenHurricane {
    namespace cuChem {
        /**
         * \brief The class of reaction table.
         */
        class reactionTable {
        public:
            enum reactionTypes : short {
                irreversibleReaction = 0,
                reversibleReaction = 1,
                nonEquilibriumReversibleReaction = 2
            };

            enum thirdBodyTypes : short {
                noThirdBody = 0,         /**< Without third-body. */
                thirdBody = 1,           /**< Only third-body. */
                UnimolecularFallOff = 2, /**< Unimolecular/recombination fall-off reactions. */
                bimolecularFallOff = 3 /**< Chemically activated bimolecular fall-off reactions. */
            };

        private:
            cu_ushort nrc_;
            cu_ushort nsp_;

        public:
            inline cu_dual cu_ushort nsp() const noexcept;
            inline cu_dual cu_ushort nrc() const noexcept;

        private:
            /** \brief Forward reaction coefficients. */
            reactionCoeffs fc_;

            /** \brief Reverse reaction coefficients. */
            reactionCoeffs rc_;

            /*!\brief The summation of stoichiometry coeffs with size of nrc.*/
            cu1DArray<cu_real> sumStoi_;

            /**
             * \brief Reaction type:
             *   0 - irreversible reaction.
             *   1 - reversible reaction.
             *   2 - non-equilibrium reversible reaction.
             */
            cu1DArray<cu_ushort> reactionType_;

            /**
             * \brief Third-body type:
             *   0 - Without third-body.
             *   1 - Only third-body.
             *   2 - Unimolecular/recombination fall-off reactions.
             *   3 - Chemically activated bimolecular fall-off reactions.
             */
            cu1DArray<cu_ushort> thirdBodyType_;

            /** \brief Forward reaction rate coefficients. */
            reactionRateCoeffs kf_;

            /** \brief Reverse reaction rate coefficients. */
            reactionRateCoeffs kr_;

            /** \brief The third-body coefficients. */
            thirdBodyCoefficients tbc_;

            /** \brief The pressure-dependent reaction coefficients. */
            pressureDependentCoefficients pdc_;

            /** \brief The table of species. */
            const speciesTable species_;

        public:
            /** \brief The table of species. */
            inline cu_dual const speciesTable &species() const;

        public:
            /**
             * \brief Null constructor.
             */
            inline cu_host reactionTable();

            /**
             * \brief Construct from interface.
             */
            inline cu_host reactionTable(const cuChemInterface &reacTableInt,
                                                const speciesTable spt);

            /**
             * \brief Construct from interface.
             */
            inline cu_host reactionTable(const cuChemInterface &reacTableInt,
                                                const speciesTable spt, const cudaStreams &streams);

            /**
             * \brief Copy constructor.
             */
            inline cu_dual reactionTable(const reactionTable &);

            /**
             * \brief Destructor.
             */
            cu_dual inline ~reactionTable() noexcept;

            /**
             * \brief Destroy arrays stored in GPU gloabl memory.
             */
            inline cu_host void destroy();

            /*!\brief Equilibrium constant i.t.o. partial pressures*/
            cu_device cu_real Kp(const cu_ushort irc, const cu_real T,
                                          const cu_real *__restrict__ GdRT) const;

            /*!
             * \brief Equilibrium constant i.t.o. molar concentration.
             * \param[in] irc - The index of the reaction
             * \param[in] T - Temperature [K]
             * \return Equilibrium constant i.t.o. molar concentration
             * \retval A real value
             */
            inline cu_device cu_real Kc(const cu_ushort irc, const cu_real T,
                                                 const cu_real *__restrict__ nGdRT) const;

            /**
             * \brief Forward reaction rate constant.
             * \param[in] irc - The index of the reaction
             * \param[in] T - Temperature [K]
             * \return Forward reaction rate constant
             * \retval A real value
             */
            inline cu_device cu_real kf(const cu_ushort irc, const cu_real T) const;

            /**
             * \brief Reverse reaction rate constant for non-equilibrium reversible reactions.
             * \param[in] irc - The index of the reaction
             * \param[in] T - Temperature [K]
             * \return Reverse reaction rate constant
             * \retval A real value
             */
            inline cu_device cu_real kr(const cu_ushort irc, const cu_real T) const;

            /**
             * \brief Reverse reaction rate constant for reversible reactions.
             * \param[in] kf - Forward reaction rate constant
             * \param[in] kc - Equilibrium constant i.t.o. molar concentration
             * \return Reverse reaction rate constant
             * \retval A real value
             */
            inline cu_device cu_real kr(const cu_real kf, const cu_real kc) const;

            /**
             * \brief Return The forward and reverse rate constant of irc_th reaction.
             * \param[in] irc - the index of reactions.
             * \param[in] T - static temperature [K].
             * \param[in] ci - molar concentrations of species [kmol/m^3].
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] kff - forward rate constant of irc_th reaction.
             * \param[out] krr - reverse rate constant of irc_th reaction.
             */
            cu_device void kfAndKr(const cu_ushort irc, const cu_real T,
                                          const cu_real *__restrict__ ci,
                                          const cu_real *__restrict__ GdRT, cu_real &kff,
                                          cu_real &krr) const;

            /*!\brief Return the net rate of irc_th reaction.
             * \param[in] irc - the index of reactions.
             * \param[in] T - static temperature [K].
             * \param[in] ci - molar concentrations of species [kmol/m^3].
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] kff - forward rate constant of irc_th reaction.
             * \param[out] krr - reverse rate constant of irc_th reaction.
             * \param[out] qfj - forward rate of irc_th reaction.
             * \param[out] qrj - reverse rate of irc_th reaction.
             * \return The net rate of irc_th reaction.
             */
            cu_device cu_real qfr(const cu_ushort irc, const cu_real T,
                                           const cu_real *__restrict__ ci,
                                           const cu_real *__restrict__ GdRT, cu_real &kff,
                                           cu_real &krr, cu_real &qf, cu_real &qr) const;

            /*!\brief Return the net rate of irc_th reaction.
             * \param[in] irc - the index of reactions.
             * \param[in] T - static temperature [K].
             * \param[in] ci - molar concentrations of species [kmol/m^3].
             * \param[in] GdRT - Gibb's free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] kff - forward rate constant of irc_th reaction.
             * \param[out] krr - reverse rate constant of irc_th reaction.
             * \param[out] qfj - forward rate of irc_th reaction.
             * \param[out] qrj - reverse rate of irc_th reaction.
             * \return The net rate of irc_th reaction.
             */
            /*cu_device cu_real qfr(const cu_ushort irc, const cu_real T,
                                           const cu_real *__restrict__ ci,
                                           const cu_real *__restrict__ GdRT, cu_real &kff,
                                           cu_real &krr, cu_real &qf, cu_real &qr,
                                           cu_real *__restrict__ dOmegaidCi) const;*/

            /*!\brief Return the net rate of irc_th reaction.
             * \param[in] irc - the index of reactions.
             * \param[in] T - static temperature [K].
             * \param[in] ci - molar concentrations of species [kmol/m^3].
             * \param[in] GdRT - Gibb's free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] kff - forward rate constant of irc_th reaction.
             * \param[out] krr - reverse rate constant of irc_th reaction.
             * \param[out] qfj - forward rate of irc_th reaction.
             * \param[out] qrj - reverse rate of irc_th reaction.
             * \return The net rate of irc_th reaction.
             */
            template<class dOmegaType>
            cu_device cu_real qfr(const cu_ushort irc, const cu_real T,
                                           const cu_real *__restrict__ ci,
                                           const cu_real *__restrict__ GdRT, cu_real &kff,
                                           cu_real &krr, cu_real &qf, cu_real &qr,
                                           dOmegaType *__restrict__ dOmegaidCi) const;

            /*!\brief Return the net rate of irc_th reaction.
             * \param[in] irc - the index of reactions.
             * \param[in] T - static temperature [K].
             * \param[in] ci - molar concentrations of species [kmol/m^3].
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] kff - forward rate constant of irc_th reaction.
             * \param[out] krr - reverse rate constant of irc_th reaction.
             * \param[out] qfj - forward rate of irc_th reaction.
             * \param[out] qrj - reverse rate of irc_th reaction.
             * \return The net rate of irc_th reaction.
             */
            cu_device cu_real qfrAllSpecies(const cu_ushort irc, const cu_real T,
                                                     const cu_real *__restrict__ ci,
                                                     const cu_real *__restrict__ GdRT,
                                                     cu_real &kff, cu_real &krr, cu_real &qf,
                                                     cu_real &qr,
                                                     cu_real *__restrict__ dOmegaidCi) const;

            /*!\brief Return the net rate of irc_th reaction.
             * \param[in] irc - the index of reactions.
             * \param[in] T - static temperature [K].
             * \param[in] ci - molar concentrations of species [kmol/m^3].
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] kff - forward rate constant of irc_th reaction.
             * \param[out] krr - reverse rate constant of irc_th reaction.
             * \return The net rate of irc_th reaction.
             */
            cu_device void dWdciAllSpecies(const cu_ushort irc, const cu_real T,
                                                  const cu_real *__restrict__ ci,
                                                  const cu_real *__restrict__ GdRT,
                                                  cu_real &kff, cu_real &krr,
                                                  cu_real *__restrict__ dOmegaidCi) const;

            /*!\brief Return the net rate of irc_th reaction.
             * \param[in] irc - the index of reactions.
             * \param[in] T - static temperature [K].
             * \param[in] ci - molar concentrations of species [kmol/m^3].
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] kff - forward rate constant of irc_th reaction.
             * \param[out] krr - reverse rate constant of irc_th reaction.
             * \return The net rate of irc_th reaction.
             */
            inline cu_device cu_real qfr(const cu_ushort irc, const cu_real T,
                                                  const cu_real *__restrict__ ci,
                                                  const cu_real *__restrict__ GdRT,
                                                  cu_real &kff, cu_real &krr) const;

            /*!\brief Return the net rate of irc_th reaction.
             * \param[in] irc - the index of reactions.
             * \param[in] T - static temperature [K].
             * \param[in] ci - molar concentrations of species [kmol/m^3].
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \return (qf-qr) - the net rate of irc_th reaction.
             */
            inline cu_device cu_real qfr(const cu_ushort irc, const cu_real T,
                                                  const cu_real *__restrict__ ci,
                                                  const cu_real *__restrict__ GdRT) const;

            /*!\brief Chemical source terms for coupled NS equations.
             * \param[in] irc - the index of reactions.
             * \param[in] T - temperature [K]
             * \param[in] c - molar concentration [kmol/m^3]
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] omegai - the chemical source term with an array size of Nsp - 1 [in molar unit]
             */
            cu_device void omegaCoupled(const cu_ushort irc, const cu_real T,
                                               const cu_real *__restrict__ ci,
                                               const cu_real *__restrict__ GdRT,
                                               cu_real *__restrict__ omegai) const;

            /*!\brief Chemical source terms for coupled NS equations.
             * \param[in] irc - the index of reactions.
             * \param[in] T - temperature [K]
             * \param[in] c - molar concentration [kmol/m^3]
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] omegai - the chemical source term with an array size of Nsp [in molar unit]
             */
            cu_device void omega(const cu_ushort irc, const cu_real T,
                                        const cu_real *__restrict__ ci,
                                        const cu_real *__restrict__ GdRT,
                                        cu_real *__restrict__ omegai) const;

            /*!\brief Chemical source terms for coupled NS equations.
             * \param[in] irc - the index of reactions.
             * \param[in] T - temperature [K]
             * \param[in] c - molar concentration [kmol/m^3]
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] omegai - the chemical source term with an array size of Nsp - 1 [in molar unit]
             */
            cu_device void omegaCoupled(const cu_ushort irc, const cu_real T,
                                               const cu_real *__restrict__ ci,
                                               const cu_real *__restrict__ GdRT,
                                               cu_real *__restrict__ omegai,
                                               cu_real *__restrict__ dOmegaidCi) const;

            /*!\brief Chemical source terms for coupled NS equations.
             * \param[in] irc - the index of reactions.
             * \param[in] T - temperature [K]
             * \param[in] c - molar concentration [kmol/m^3]
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] omegai - the chemical source term with an array size of Nsp - 1 [in molar unit]
             */
            cu_device void omegaCoupled2(const cu_ushort irc, const cu_real T,
                                                const cu_real *__restrict__ ci,
                                                const cu_real *__restrict__ GdRT,
                                                cu_real *__restrict__ omegai,
                                                cu_real *__restrict__ dOmegaidCi) const;

            /*!\brief Chemical source terms for coupled NS equations.
             * \param[in] irc - the index of reactions.
             * \param[in] T - temperature [K]
             * \param[in] c - molar concentration [kmol/m^3]
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] omegai - the chemical source term with an array size of Nsp - 1 [in molar unit]
             */
            /*cu_device void omegaCoupled3(const cu_ushort irc, const cu_real T,
                                                const cu_real *__restrict__ ci,
                                                const cu_real *__restrict__ GdRT,
                                                cu_real *__restrict__ omegai,
                                                cu_real *__restrict__ dOmegaidCi) const;*/

            /*!\brief Chemical source terms for coupled NS equations.
             * \param[in] irc - the index of reactions.
             * \param[in] T - temperature [K]
             * \param[in] c - molar concentration [kmol/m^3]
             * \param[in] GdRT - Gibb¡¯s free energy: \f$ - \frac{{G_i^0}}{{{R_0}T}} = \frac{{TS_i^0 - H_i^0}}{{{R_0}T}}\f$.
             * \param[out] omegai - the chemical source term with an array size of Nsp - 1 [in molar unit]
             */
            template <class dOmegaType>
            cu_device void omegaCoupled3(const cu_ushort irc, const cu_real T,
                                                const cu_real *__restrict__ ci,
                                                const cu_real *__restrict__ GdRT,
                                                cu_real *__restrict__ omegai,
                                                dOmegaType *__restrict__ dOmegaidCi) const;

            /*!\brief Chemical source terms for coupled NS equations.
             * \param[in] isp - the index of species.
             * \param[in] omegai - the chemical source term with an array size of Nsp - 1 [in molar unit]
             * \param[in] rhoRef - density reference value
             * \param[in] timeRef - time reference value
             * \param[out] Rii - The chemical source terms of species
             */
            inline cu_device void Ri(const cu_ushort isp,
                                            const cu_real *__restrict__ omegai,
                                            cu_real *__restrict__ Rii) const;
        };
    } // namespace cuChem
} // namespace OpenHurricane

#include "reactionTableCUDA.inl"
#endif // CUDA_PARALLEL
