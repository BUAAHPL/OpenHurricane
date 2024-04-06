/*!
 * \file equationOfState.hpp
 * \brief Header of base class of equation of state.
 *       The subroutines and functions are in the <i>equationOfState.cpp</i> file.
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

#include "dataStructure.hpp"
#include "objectFactory.hpp"
#include "speciesList.hpp"

namespace OpenHurricane {
    /*!\brief The base class of equation of state.*/
    class equationOfState {
    private:
        /*!\brief Reference to the species list.*/
        const speciesList &species_;

    public:
        declareClassNames;
        declareObjFty(equationOfState, controller,
                      (const speciesList &sp, const controller &eosCont), (sp, eosCont));

        /*!\brief Construct from components.*/
        inline equationOfState(const speciesList &_species, const controller &cont)
            : species_(_species) {}

        /*!\brief Return a clone.*/
        virtual hur_nodiscard uniquePtr<equationOfState> clone() const = 0;

        /*!\brief Select null constructed.*/
        hur_nodiscard static uniquePtr<equationOfState> creator(const speciesList &_species,
                                                  const controller &cont);

        /*!\brief Destructor.*/
        virtual ~equationOfState() noexcept {}

        /*!\brief Return the species list.*/
        hur_nodiscard inline const speciesList &species() const noexcept { return species_; }

        hur_nodiscard inline integer size() const noexcept { return species_.size(); }

        /**\brief Is the equation of state incompressible, i.e. rho != f(p).*/
        hur_nodiscard virtual bool isIncompressible() const noexcept = 0;

        /**\brief Is the equation of state isochoric, i.e. rho = const.*/
        hur_nodiscard virtual bool isIsochoric() const noexcept = 0;

        /** \brief Return the density of specie i. [kg/m^3]*/
        hur_nodiscard virtual real rhoi(const real pi, const real Ti, const integer i) const = 0;

        /** \brief Return the density of mixture. [kg/m^3]*/
        hur_nodiscard virtual real rhom(const real pm, const real Tm,
                                        const realArray &_yi) const = 0;

        /** \brief Return the density of mixture. [kg/m^3]*/
        hur_nodiscard virtual real rhom(const real pm, const real Tm,
                                        const PtrList<cellRealArray> &_yi,
                                        const integer celli) const = 0;

        /** \brief Return the pressure of specie i. [Pa]*/
        hur_nodiscard virtual real pi(const real rhoi, const real Ti, const integer i) const = 0;

        /** \brief Return the pressure of mixture. [Pa]*/
        hur_nodiscard virtual real pm(const real rhom, const real Tm,
                                      const realArray &_yi) const = 0;

        /** \brief Return the pressure of mixture. [Pa]*/
        hur_nodiscard virtual real pm(const real rhom, const real Tm,
                                      const PtrList<cellRealArray> &_yi,
                                      const integer celli) const = 0;

        /** \brief Return the temperature of specie i. [K]*/
        hur_nodiscard virtual real Ti(const real rhoi, const real pi, const integer i) const = 0;

        /** \brief Return the temperature of mixture. [Pa]*/
        hur_nodiscard virtual real Tm(const real rhom, const real pm,
                                      const realArray &_yi) const = 0;

        /** \brief Return the temperature of mixture. [Pa]*/
        hur_nodiscard virtual real Tm(const real rhom, const real pm,
                                      const PtrList<cellRealArray> &_yi,
                                      const integer celli) const = 0;

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg] ��ʣ���ʣ�
         *         by given pressure
         * \param[in] pi - pressure for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         */
        hur_nodiscard virtual real hi_p(const real p, const real T, const integer i) const = 0;

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg] ��ʣ���ʣ�
         *         by given pressure
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         */
        hur_nodiscard virtual real hm_p(const real p, const real T, const realArray &yi) const = 0;

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg] ��ʣ���ʣ�
         *         by given pressure
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         */
        hur_nodiscard virtual real hm_p(const real p, const real T,
                                        const PtrList<cellRealArray> &_yi,
                                        const integer celli) const = 0;

        /** \brief Return molar departure enthalpy of specie i. [J/kmol] ��ʣ���ʣ�
         *         by given pressure
         * \param[in] pi - pressure for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         */
        hur_nodiscard virtual real Hi_p(const real p, const real T, const integer i) const = 0;

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg] ��ʣ���ʣ�
         *         by given rho (density)
         * \param[in] rhoi - density for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         */
        hur_nodiscard virtual real hi_rho(const real rho, const real T, const integer i) const = 0;

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg] ��ʣ���ʣ�
         *         by given rho (density)
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         */
        hur_nodiscard virtual real hm_rho(const real rho, const real T,
                                          const realArray &yi) const = 0;

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg] ��ʣ���ʣ�
         *         by given rho (density)
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         */
        hur_nodiscard virtual real hm_rho(const real rho, const real T,
                                          const PtrList<cellRealArray> &_yi,
                                          const integer celli) const = 0;

        /** \brief Return molar departure enthalpy of specie i. [J/kmol] ��ʣ���ʣ�
         *         by given rho (density)
         * \param[in] rhoi - density for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         */
        hur_nodiscard virtual real Hi_rho(const real rho, const real T, const integer i) const = 0;

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given pressure
         * \param[in] pi - pressure for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         */
        hur_nodiscard virtual real cpi_p(const real p, const real T, const integer i) const = 0;

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given pressure
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         */
        hur_nodiscard virtual real cpm_p(const real p, const real T, const realArray &yi) const = 0;

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given pressure
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         */
        hur_nodiscard virtual real cpm_p(const real p, const real T,
                                         const PtrList<cellRealArray> &_yi,
                                         const integer celli) const = 0;

        /** \brief Return molar departure specific heat at constant pressure [J/(kmol K)]
         *         by given pressure
         * \param[in] pi - pressure for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         */
        hur_nodiscard virtual real Cpi_p(const real p, const real T, const integer i) const = 0;

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kmol K)]
         *         by given pressure
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         */
        hur_nodiscard virtual real Cpm_p(const real p, const real T, const realArray &xi) const = 0;

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given rho (density).
         * \param[in] rhoi - density for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         */
        hur_nodiscard virtual real cpi_rho(const real rho, const real T, const integer i) const = 0;

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given rho (density).
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species fractions field.
         */
        hur_nodiscard virtual real cpm_rho(const real rho, const real T,
                                           const realArray &yi) const = 0;

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given rho (density).
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species fractions field.
         */
        hur_nodiscard virtual real cpm_rho(const real rho, const real T,
                                           const PtrList<cellRealArray> &_yi,
                                           const integer celli) const = 0;

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kmol K)]
         *         by given rho (density).
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] xi - species molar fractions field.
         */
        hur_nodiscard virtual real Cpm_rho(const real rho, const real T,
                                           const realArray &yi) const = 0;

        /** \brief Return molar departure specific heat at constant pressure [J/(kmol K)]
         *         by given rho (density).
         * \param[in] rhoi - density for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         */
        hur_nodiscard virtual real Cpi_rho(const real rho, const real T, const integer i) const = 0;

        /**
         * \brief Return specific entropy [J/(kg K)].
         *        by given pressure.
         * \param[in] pi - pressure for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard virtual real si_p(const real p, const real T, const integer i) const = 0;

        /**
         * \brief Return specific entropy [J/(kg K)].
         *        by given pressure.
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fraction field.
         *
         */
        hur_nodiscard virtual real sm_p(const real p, const real T, const realArray &yi) const = 0;

        /**
         * \brief Return specific entropy [J/(kg K)].
         *        by given pressure.
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fraction field.
         *
         */
        hur_nodiscard virtual real sm_p(const real p, const real T,
                                        const PtrList<cellRealArray> &_yi,
                                        const integer celli) const = 0;

        /**
         * \brief Return specific entropy [J/(kg K)].
         *        by given density.
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard virtual real si_rho(const real rho, const real T, const integer i) const = 0;

        /**
         * \brief Return specific entropy [J/(kg K)].
         *        by given density.
         * \param[in] rho - density for mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         *
         */
        hur_nodiscard virtual real sm_rho(const real rho, const real T,
                                          const realArray &Yi) const = 0;

        /**
         * \brief Return specific entropy [J/(kg K)].
         *        by given density.
         * \param[in] rho - density for mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         *
         */
        hur_nodiscard virtual real sm_rho(const real rho, const real T,
                                          const PtrList<cellRealArray> &_yi,
                                          const integer celli) const = 0;

        /**
         * \brief Return molar entropy [J/(kmol K)].
         *        by given pressure.
         * \param[in] pi - pressure for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard virtual real Si_p(const real p, const real T, const integer i) const = 0;

        /**
         * \brief Return molar entropy [J/(kmol K)].
         *        by given pressure.
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] xi - species molar fraction field.
         *
         */
        hur_nodiscard virtual real Sm_p(const real p, const real T, const realArray &xi) const = 0;

        /**
         * \brief Return molar entropy [J/(kmol K)].
         *        by given density.
         * \param[in] rhoi - density for mixtures.
         * \param[in] Ti - temperature for mixtures.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard virtual real Si_rho(const real rhoi, const real Ti,
                                          const integer i) const = 0;

        /**
         * \brief Return molar entropy [J/(kmol K)].
         *        by given density.
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] xi - species molar fractions field.
         *
         */
        hur_nodiscard virtual real Sm_rho(const real rho, const real T,
                                          const realArray &xi) const = 0;

        /**
         * \brief Return compression factor: Z = p/(rho * R * T).
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard virtual real Zi(const real pi, const real Ti, const integer i) const = 0;

        /**
         * \brief Return compression factor: Z = p/(rho * R * T).
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard virtual real Zm(const real pm, const real Tm, const realArray &_yi) const = 0;

        /**
         * \brief Return compression factor: Z = p/(rho * R * T).
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard virtual real Zm(const real pm, const real Tm,
                                      const PtrList<cellRealArray> &_yi,
                                      const integer celli) const = 0;

        /**\brief Return compressibility: rho/p [s^2/m^2].
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard virtual real psii(const real pi, const real Ti, const integer i) const = 0;

        /**\brief Return compressibility: rho/p [s^2/m^2].
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard virtual real psim(const real pm, const real Tm,
                                        const realArray &_yi) const = 0;

        /**\brief Return compressibility: rho/p [s^2/m^2].
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard virtual real psim(const real pm, const real Tm,
                                        const PtrList<cellRealArray> &_yi,
                                        const integer celli) const = 0;

        /**
         * \brief Return (cpi - cvi) [J/(kg K)].
         *        by given pressure.
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard virtual real cpMcvi_p(const real pi, const real Ti,
                                            const integer i) const = 0;

        /**
         * \brief Return (cpi - cvi) [J/(kg K)].
         *        by given pressure.
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard virtual real cpMcvm_p(const real pm, const real Tm,
                                            const realArray &_yi) const = 0;

        /**
         * \brief Return (cpi - cvi) [J/(kg K)].
         *        by given pressure.
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard virtual real cpMcvm_p(const real pm, const real Tm,
                                            const PtrList<cellRealArray> &_yi,
                                            const integer celli) const = 0;

        /**
         * \brief Return (cpi - cvi) [J/(kg K)].
         *        by given density.
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard virtual real cpMcvi_rho(const real rhoi, const real Ti,
                                              const integer i) const = 0;

        /**
         * \brief Return (Cpi - Cvi) [J/(kmol K)].
         *        by given pressure.
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard virtual real CpMCvi_p(const real pi, const real Ti,
                                            const integer i) const = 0;

        /**
         * \brief Return (Cpi - Cvi) [J/(kmol K)].
         *        by given density.
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard virtual real CpMCvi_rho(const real rhoi, const real Ti,
                                              const integer i) const = 0;

        hur_nodiscard virtual real DRhoDp(const real rhoi, const real Ti,
                                          const integer i) const = 0;
        hur_nodiscard virtual real DRhoDT(const real rhoi, const real Ti,
                                          const integer i) const = 0;
        hur_nodiscard virtual real DpDT(const real rhoi, const real Ti, const integer i) const = 0;
        hur_nodiscard virtual real DpDRho(const real rhoi, const real Ti,
                                          const integer i) const = 0;

        hur_nodiscard virtual real DRhoDp(const real p, const real rho, const real T,
                                          const realArray &y) const = 0;

        hur_nodiscard virtual real DRhoDp(const real p, const real rho, const real T,
                                          const PtrList<cellRealArray> &_yi,
                                          const integer celli) const = 0;

        hur_nodiscard virtual real DRhoDT(const real p, const real rho, const real T,

                                          const realArray &y) const = 0;

        hur_nodiscard virtual real DRhoDT(const real p, const real rho, const real T,
                                          const PtrList<cellRealArray> &_yi,
                                          const integer celli) const = 0;

        hur_nodiscard virtual real DRhoInverseDT(const real p, const real rho, const real T,

                                                 const realArray &y) const = 0;

        hur_nodiscard virtual real DpDT(const real p, const real rho, const real T,
                                        const realArray &y) const = 0;

        hur_nodiscard virtual real DpDT(const real p, const real rho, const real T,
                                        const PtrList<cellRealArray> &_yi,
                                        const integer celli) const = 0;

        hur_nodiscard virtual real DpDRho(const real p, const real rho, const real T,
                                          const realArray &y) const = 0;

        hur_nodiscard virtual real DpDRho(const real p, const real rho, const real T,
                                          const PtrList<cellRealArray> &_yi,
                                          const integer celli) const = 0;

        hur_nodiscard virtual real gammaThCorrect(const real rhoi, const real Ti,
                                                  const integer i) const {
            return real(1.0);
        }

        hur_nodiscard virtual real gammaThCorrectM(const real p, const real rho, const real T,
                                                   const realArray &yi) const {
            return real(1.0);
        }

        hur_nodiscard virtual real gammaThCorrectM(const real p, const real rho, const real T,
                                                   const PtrList<cellRealArray> &_yi,
                                                   const integer celli) const {
            return real(1.0);
        }
    };
} // namespace OpenHurricane
