/*!
 * \file perfectGas.hpp
 * \brief Header of equation of state of perfect gas.
 *       The subroutines and functions are in the <i>perfectGas.cpp</i> file.
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

#include "commonInclude.hpp"
#include "equationOfState.hpp"
#include "real.hpp"
#include "smartPointer.hpp"
#include "string.hpp"

namespace OpenHurricane {
    /**
     * \brief The class of the equation of state of perfect gas.
     *       p = rho*R*T
     */
    class perfectGas : public equationOfState {
    public:
        /**\brief Is the equation of state incompressible, i.e. rho != f(p).*/
        static constexpr bool incompressible = false;

        /**\brief Is the equation of state isochoric, i.e. rho = const.*/
        static constexpr bool isochoric = false;

    public:
        /** \brief Return the bcType name*/
        declareClassNames;

        /** \brief Construct from components*/
        inline perfectGas(const speciesList &_species, const controller &cont)
            : equationOfState(_species, cont) {}

        /** \brief Return a clone*/
        inline virtual hur_nodiscard uniquePtr<equationOfState> clone() const {
            return uniquePtr<equationOfState>(new perfectGas(*this));
        }

        /*!\brief Destructor.*/
        virtual ~perfectGas() noexcept {}

        /**\brief Is the equation of state incompressible, i.e. rho != f(p).*/
        hur_nodiscard virtual bool isIncompressible() const noexcept {
            return perfectGas::incompressible;
        }

        /**\brief Is the equation of state isochoric, i.e. rho = const.*/
        hur_nodiscard virtual bool isIsochoric() const noexcept { return perfectGas::isochoric; }

        /** \brief Return the density of specie i [kg/m^3]*/
        hur_nodiscard inline real rhoi(const real pi, const real Ti,
                                       const integer i) const noexcept {
            return pi / (species()[i].Ri() * Ti);
        }

        /** \brief Return the density of mixture. [kg/m^3]*/
        hur_nodiscard virtual real rhom(const real pm, const real Tm,
                                        const realArray &_yi) const noexcept;

        /** \brief Return the density of mixture. [kg/m^3]*/
        hur_nodiscard virtual real rhom(const real pm, const real Tm,
                                        const PtrList<cellRealArray> &_yi,
                                        const integer celli) const noexcept;

        /** \brief Return the pressure of specie i. [Pa]*/
        hur_nodiscard inline virtual real pi(const real rhoi, const real Ti,
                                             const integer i) const noexcept {
            return rhoi * (species()[i].Ri() * Ti);
        }

        /** \brief Return the pressure of mixture. [Pa]*/
        hur_nodiscard virtual real pm(const real rhom, const real Tm,
                                      const realArray &_yi) const noexcept;

        /** \brief Return the pressure of mixture. [Pa]*/
        hur_nodiscard virtual real pm(const real rhom, const real Tm,
                                      const PtrList<cellRealArray> &_yi,
                                      const integer celli) const noexcept;

        /** \brief Return the temperature of specie i. [K]*/
        hur_nodiscard inline virtual real Ti(const real rhoi, const real pi,
                                             const integer i) const noexcept {
            return pi / (rhoi * species().Ri(i));
        }

        /** \brief Return the temperature of mixture. [Pa]*/
        hur_nodiscard inline virtual real Tm(const real rhom, const real pm,
                                             const realArray &_yi) const noexcept {
            real Rm = species().Rm(_yi);
            return pm / (rhom * Rm);
        }

        /** \brief Return the temperature of mixture. [Pa]*/
        hur_nodiscard inline virtual real Tm(const real rhom, const real pm,
                                             const PtrList<cellRealArray> &_yi,
                                             const integer celli) const noexcept {
            real Rm = species().Rm(_yi, celli);
            return pm / (rhom * Rm);
        }

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg]
         *         by given pressure
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard inline virtual real hi_p(const real pi, const real Ti,
                                               const integer i) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg]
         *         by given pressure
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         */
        hur_nodiscard inline virtual real hm_p(const real p, const real T,
                                               const realArray &yi) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg]
         *         by given pressure
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         */
        hur_nodiscard inline virtual real hm_p(const real p, const real T,
                                               const PtrList<cellRealArray> &_yi,
                                               const integer celli) const noexcept {
            return real(0.0);
        }

        /** \brief Return molar departure enthalpy of specie i. [J/kmol]
         *         by given pressure
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard inline virtual real Hi_p(const real pi, const real Ti,
                                               const integer i) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg]
         *         by given rho (density)
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard inline virtual real hi_rho(const real rhoi, const real Ti,
                                                 const integer i) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg]
         *         by given rho (density)
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         */
        hur_nodiscard virtual real hm_rho(const real rho, const real T,
                                          const realArray &yi) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure enthalpy of specie i. [J/kg]
         *         by given rho (density)
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         */
        hur_nodiscard virtual real hm_rho(const real rho, const real T,
                                          const PtrList<cellRealArray> &_yi,
                                          const integer celli) const noexcept {
            return real(0.0);
        }

        /** \brief Return molar departure enthalpy of specie i. [J/kmol]
         *         by given rho (density)
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard inline virtual real Hi_rho(const real rhoi, const real Ti,
                                                 const integer i) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given pressure
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard inline virtual real cpi_p(const real p, const real Ti,
                                                const integer i) const noexcept {
            return real(0.0);
        }

        /** \brief Return molar departure specific heat at constant pressure [J/(kmol K)]
         *         by given pressure
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard inline virtual real Cpi_p(const real p, const real Ti,
                                                const integer i) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given pressure
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         */
        hur_nodiscard virtual real cpm_p(const real p, const real T,
                                         const realArray &yi) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given pressure
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         */
        hur_nodiscard virtual real cpm_p(const real p, const real T,
                                         const PtrList<cellRealArray> &_yi,
                                         const integer celli) const noexcept {
            return real(0.0);
        }

        hur_nodiscard virtual real Cpm_p(const real p, const real T,
                                         const realArray &xi) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given rho (density).
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard inline virtual real cpi_rho(const real rho, const real Ti,
                                                  const integer i) const noexcept {
            return real(0.0);
        }

        /** \brief Return molar departure specific heat at constant pressure [J/(kmol K)]
         *         by given rho (density).
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard inline virtual real Cpi_rho(const real rho, const real Ti,
                                                  const integer i) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given rho (density).
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species fractions field.
         */
        hur_nodiscard virtual real cpm_rho(const real rho, const real T,
                                           const realArray &yi) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kmol K)]
         *         by given rho (density).
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] xi - species molar fractions field.
         */
        hur_nodiscard virtual real Cpm_rho(const real rho, const real T,
                                           const realArray &xi) const noexcept {
            return real(0.0);
        }

        /** \brief Return specific (mass-based) departure specific heat at constant pressure [J/(kg K)]
         *         by given rho (density).
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species fractions field.
         */
        hur_nodiscard virtual real cpm_rho(const real rho, const real T,
                                           const PtrList<cellRealArray> &_yi,
                                           const integer celli) const noexcept {
            return real(0.0);
        }

        /**
         * \brief Return specific entropy [J/(kg K)].
         *        by given pressure.
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard inline virtual real si_p(const real pi, const real Ti,
                                               const integer i) const {
            return -species()[i].Ri() * log(max(pi, veryTiny) / constant::physicalConstant::Patm);
        }

        /**
         * \brief Return specific entropy [J/(kg K)].
         *        by given pressure.
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fraction field.
         *
         */
        hur_nodiscard virtual real sm_p(const real p, const real T, const realArray &yi) const;

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
                                        const integer celli) const;

        /**
         * \brief Return specific entropy [J/(kg K)].
         *        by given density.
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard inline virtual real si_rho(const real rhoi, const real Ti,
                                                 const integer i) const {
            real pk = this->pi(rhoi, Ti, i);
            return si_p(pk, Ti, i);
        }

        /**
         * \brief Return specific entropy [J/(kg K)].
         *        by given density.
         * \param[in] rho - density for species mixturei.
         * \param[in] T - temperature for species mixture.
         * \param[in] yi - species mass fractions field.
         *
         */
        hur_nodiscard virtual real sm_rho(const real rho, const real T, const realArray &Yi) const;

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
                                          const integer celli) const;

        /**
         * \brief Return molar entropy [J/(kmol K)].
         *        by given pressure.
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard inline virtual real Si_p(const real pi, const real Ti,
                                               const integer i) const {
            return -constant::physicalConstant::Ru *
                   log(max(pi, veryTiny) / constant::physicalConstant::Patm);
        }

        /**
         * \brief Return molar entropy [J/(kmol K)].
         *        by given pressure.
         * \param[in] p - pressure for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] xi - species molar fraction field.
         *
         */
        hur_nodiscard virtual real Sm_p(const real p, const real T, const realArray &xi) const;

        /**
         * \brief Return molar entropy [J/(kmol K)].
         *        by given density.
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard inline virtual real Si_rho(const real rhoi, const real Ti,
                                                 const integer i) const {
            real pk = this->pi(rhoi, Ti, i);
            return Si_p(pk, Ti, i);
        }

        /**
         * \brief Return molar entropy [J/(kmol K)].
         *        by given density.
         * \param[in] rho - density for species mixture.
         * \param[in] T - temperature for species mixture.
         * \param[in] xi - species molar fractions field.
         *
         */
        hur_nodiscard virtual real Sm_rho(const real rho, const real T, const realArray &xi) const;

        /**
         * \brief Return compression factor: Z = p/(rho * R * T).
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard inline virtual real Zi(const real pi, const real Ti,
                                             const integer i) const noexcept {
            return real(1.0);
        }

        /**
         * \brief Return compression factor: Z = p/(rho * R * T).
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard inline virtual real Zm(const real pm, const real Tm,
                                             const realArray &_yi) const noexcept {
            return real(1.0);
        }

        /**
         * \brief Return compression factor: Z = p/(rho * R * T).
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard virtual real Zm(const real pm, const real Tm,
                                      const PtrList<cellRealArray> &_yi,
                                      const integer celli) const noexcept {
            return real(1.0);
        }

        /**\brief Return compressibility: rho/p [s^2/m^2].
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         */
        hur_nodiscard inline virtual real psii(const real pi, const real Ti,
                                               const integer i) const noexcept {
            return real(1.0) / (species()[i].Ri() * Ti);
        }

        /**\brief Return compressibility: rho/p [s^2/m^2].
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard virtual real psim(const real pm, const real Tm,
                                        const realArray &_yi) const noexcept {
            return real(1.0) / (species().Rm(_yi) * Tm);
        }

        /**\brief Return compressibility: rho/p [s^2/m^2].
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard virtual real psim(const real pm, const real Tm,
                                        const PtrList<cellRealArray> &_yi,
                                        const integer celli) const noexcept {
            return real(1.0) / (species().Rm(_yi, celli) * Tm);
        }

        /**
         * \brief Return (cpi - cvi) [J/(kg K)].
         *        by given pressure.
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard inline virtual real cpMcvi_p(const real pi, const real Ti,
                                                   const integer i) const noexcept {
            return species()[i].Ri();
        }

        /**
         * \brief Return (cpi - cvi) [J/(kg K)].
         *        by given pressure.
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard inline virtual real cpMcvm_p(const real pm, const real Tm,
                                                   const realArray &_yi) const noexcept {
            return species().Rm(_yi);
        }

        /**
         * \brief Return (cpi - cvi) [J/(kg K)].
         *        by given pressure.
         * \param[in] pm - pressure for mixture.
         * \param[in] Tm - temperature for mixture.
         * \param[in] _yi - mass fraction of species i.
         */
        hur_nodiscard virtual real cpMcvm_p(const real pm, const real Tm,
                                            const PtrList<cellRealArray> &_yi,
                                            const integer celli) const noexcept {
            return species().Rm(_yi, celli);
        }

        /**
         * \brief Return (cpi - cvi) [J/(kg K)].
         *        by given density.
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard inline virtual real cpMcvi_rho(const real rhoi, const real Ti,
                                                     const integer i) const noexcept {
            return species()[i].Ri();
        }

        /**
         * \brief Return (Cpi - Cvi) [J/(kmol K)].
         *        by given pressure.
         * \param[in] pi - pressure for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard inline virtual real CpMCvi_p(const real pi, const real Ti,
                                                   const integer i) const noexcept {
            return constant::physicalConstant::Ru;
        }

        /**
         * \brief Return (Cpi - Cvi) [J/(kmol K)].
         *        by given density.
         * \param[in] rhoi - density for species i.
         * \param[in] Ti - temperature for species i.
         * \param[in] i - species i.
         *
         */
        hur_nodiscard inline virtual real CpMCvi_rho(const real rhoi, const real Ti,
                                                     const integer i) const noexcept {
            return constant::physicalConstant::Ru;
        }

        /*!\brief Acoustic velocity.*/
        hur_nodiscard inline real c(const real gam, const real p, const real rho) const noexcept {
            return sqrt(gam * p / rho);
        }

        hur_nodiscard inline virtual real DRhoDp(const real rhoi, const real Ti,
                                                 const integer i) const noexcept {
            return (species()[i].Ri() * Ti);
        }

        hur_nodiscard inline virtual real DRhoDT(const real rhoi, const real Ti,
                                                 const integer i) const noexcept {
            return -pi(rhoi, Ti, i) / (species()[i].Ri() * sqr(Ti));
        }

        hur_nodiscard inline virtual real DpDT(const real rhoi, const real Ti,
                                               const integer i) const noexcept {
            return species()[i].Ri() * rhoi;
        }

        hur_nodiscard inline virtual real DpDRho(const real rhoi, const real Ti,
                                                 const integer i) const noexcept {
            return species()[i].Ri() * Ti;
        }

        hur_nodiscard inline virtual real DRhoDp(const real p, const real rho, const real T,
                                                 const realArray &y) const noexcept {
            return species().Rm(y) * T;
        }

        hur_nodiscard virtual real DRhoDp(const real p, const real rho, const real T,
                                          const PtrList<cellRealArray> &_yi,
                                          const integer celli) const noexcept {
            return inv(species().Rm(_yi, celli) * T);
        }

        hur_nodiscard inline virtual real DRhoDT(const real p, const real rho, const real T,
                                                 const realArray &y) const noexcept {
            return -p / (species().Rm(y) * sqr(T));
        }

        hur_nodiscard virtual real DRhoDT(const real p, const real rho, const real T,
                                          const PtrList<cellRealArray> &_yi,
                                          const integer celli) const noexcept {
            return -p / (species().Rm(_yi, celli) * sqr(T));
        }

        hur_nodiscard virtual real DRhoInverseDT(const real p, const real rho, const real T,
                                                 const realArray &y) const noexcept {
            return species().Rm(y) / p;
        }

        hur_nodiscard inline virtual real DpDT(const real p, const real rho, const real T,
                                               const realArray &y) const noexcept {
            return species().Rm(y) * rho;
        }

        hur_nodiscard virtual real DpDT(const real p, const real rho, const real T,
                                        const PtrList<cellRealArray> &_yi,
                                        const integer celli) const noexcept {
            return species().Rm(_yi, celli) * rho;
        }

        hur_nodiscard inline virtual real DpDRho(const real p, const real rho, const real T,
                                                 const realArray &y) const noexcept {
            return species().Rm(y) * T;
        }

        hur_nodiscard virtual real DpDRho(const real p, const real rho, const real T,
                                          const PtrList<cellRealArray> &_yi,
                                          const integer celli) const noexcept {
            return species().Rm(_yi, celli) * T;
        }
    };

} // namespace OpenHurricane
