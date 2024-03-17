/*!
 * \file thermo.hpp
 * \brief Header of thermo properties.
 *       The subroutines and functions are in the <i>thermo.cpp</i> file.
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

#include "OpenHurricane.hpp"
#include "equationOfState.hpp"
#include "logFile.hpp"

namespace OpenHurricane {
    /*!\brief The class of thermo properties.*/
    class thermo {
    public:
        using ThermoType = thermo;

        /**\brief The phase bcType of species.*/
        enum class phaseType : short { solid, liquid, gas };

    private:

        /*!\brief Convergence tolerance of energy to temperature inversion functions.*/
        static const real tol_;

        /*!\brief Max number of iterations in energy to temperature inversion functions.*/
        static const int maxIter_;

        /*!\brief Hold reference to the equation of state.*/
        const equationOfState *eosPtr_;

        /**\brief The phase of species.*/
        phaseType phase_;

    protected:

        /*!\biref The index of the species.*/
        integer speciesIndex_;

    public:
        declareClassNames;

        declareObjFty(thermo, controller,
                      (const controller &cont, const equationOfState &EOS, const integer id),
                      (cont, EOS, id));

        thermo() = delete;

        thermo(const controller &cont, const equationOfState &EOS, const integer id);

        thermo(const equationOfState &EOS, const integer id);

        thermo(const equationOfState &EOS, const integer id, const phaseType pt);

        thermo(const equationOfState &EOS, const integer id, const std::string &pt);

        thermo(const thermo &t);
        thermo &operator=(const thermo &);

        static uniquePtr<thermo> creator(const controller &cont, const equationOfState &EOS,
                                         const integer id);

        /*!\brief Return a clone.*/
        virtual hur_nodiscard uniquePtr<thermo> clone() const = 0;

        /*!\brief Destructor.*/
        inline virtual ~thermo() noexcept;

        void setEos(const equationOfState *eos) { eosPtr_ = eos; }

        /*!
         * \brief Return the temperature corresponding to the value of the
         *  thermodynamic property f, given the function f = F(p, T)
         *  and dF(p, T)/dT.
         */
        real T(real f, real p, real T0, integer &iFlag,
               real (thermo::*F)(const real, const real) const,
               real (thermo::*dFdT)(const real, const real) const,
               real (thermo::*limit)(const real, integer &) const) const;

        const equationOfState &eos() const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] at standard pressure (1 atm).*/
        inline virtual real cp0(const real T) const = 0;

        inline virtual real inteCp0dT(const real T1, const real T2) const;

        /*!\brief The standard-state molar heat capacity at constant pressure.*/
        inline real Cp0(const real T) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)].*/
        inline virtual real cp_p(const real pi, const real T) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)].*/
        inline virtual real Dcp_pDT(const real pi, const real T) const;

        /*!\brief The molar heat capacity at constant pressure.*/
        inline real Cp_p(const real pi, const real T) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)].*/
        inline real DCp_pDT(const real pi, const real T) const;

        /*!\brief Heat capacity at constant pressure [J/(kg K)].*/
        inline virtual real cp_rho(const real rhoi, const real T) const;

        /*!\brief The molar heat capacity at constant pressure.*/
        inline real Cp_rho(const real rhoi, const real T) const;

        /*!\brief Heat capacity at constant volume [J/(kg K)] at standard pressure (1 atm).*/
        inline virtual real cv0(const real T) const;

        /*!\brief The standard-state molar heat capacity at constant volume.*/
        inline real Cv0(const real T) const;

        /*!\brief Heat capacity at constant volume [J/(kg K)].*/
        inline virtual real cv_p(const real pi, const real T) const;

        /*!\brief Heat capacity at constant volume [J/(kg K)].*/
        inline virtual real Dcv_pDT(const real pi, const real T) const;

        /*!\brief Molar heat capacity at constant volume [J/(kmol K)].*/
        inline virtual real Cv_p(const real pi, const real T) const;

        /*!\brief Molar heat capacity at constant volume [J/(kmol K)].*/
        inline virtual real DCv_pDT(const real pi, const real T) const;

        /*!\brief Heat capacity at constant volume [J/(kg K)].*/
        inline virtual real cv_rho(const real rhoi, const real T) const;

        /*!\brief Ratio of specific heat capacities.*/
        inline virtual real gamma_p(const real pi, const real T) const;

        /*!\brief Ratio of specific heat capacities.*/
        inline virtual real gamma_rho(const real rhoi, const real T) const;

        /*!\brief Thermo ratio of specific heat capacities for calculation of ascoutic speed.*/
        inline virtual real gammaThCorrect(const real rhoi, const real T) const;

        /*!\brief The standard-state internal energy.*/
        inline real u0(const real T) const;

        /*!\breif The standard-state molar internal energy.*/
        inline real U0(const real T) const;

        /*!\brief Absolute Enthalpy [J/kg] at standard pressure.*/
        inline virtual real ha0(const real T) const = 0;

        /*!\brief Absolute Enthalpy [J/kg] at standard pressure.*/
        inline virtual real Dha0DT(const real T) const = 0;

        /*!\brief The standard-state molar absolute enthalpy.*/
        inline real Ha0(const real T) const;

        /*!\brief The standard-state molar absolute enthalpy.*/
        inline real DHa0DT(const real T) const;

        /*!\brief Absolute Enthalpy [J/kg].*/
        inline virtual real ha_p(const real pi, const real T) const;

        /*!\brief Absolute Enthalpy [J/kg].*/
        inline virtual real ha_rho(const real rhoi, const real T) const;

        /*!\brief Absolute Enthalpy [J/kmol].*/
        inline virtual real Ha_p(const real pi, const real T) const;

        /*!\brief Absolute Enthalpy [J/kmol].*/
        inline virtual real Ha_rho(const real rhoi, const real T) const;

        /*!\brief Sensible enthalpy [J/kg].*/
        inline virtual real hs_p(const real pi, const real T) const;

        /*!\brief Molar sensible enthalpy [J/kmol].*/
        inline virtual real Hs_p(const real pi, const real T) const;

        /*!\brief Sensible enthalpy [J/kg].*/
        inline virtual real hs_rho(const real rhoi, const real T) const;

        /*!\brief Molar sensible enthalpy [J/kg].*/
        inline virtual real Hs_rho(const real rhoi, const real T) const;

        /*!\brief Chemical enthalpy (Enthalpy of formation) [J/kg].*/
        inline virtual real hc() const = 0;

        /*!\brief Molar chemical enthalpy (Enthalpy of formation) [J/kmol].*/
        inline real Hc() const;

        /*!\brief Absolute internal energy [J/kg].*/
        inline virtual real ea(const real rho, const real p, const real T) const;

        /*!\brief Entropy [J/(kg K)] at standard pressure.*/
        inline virtual real s0(const real T) const = 0;

        /*!\brief Molar entropy [J/(kmol K)] at standard pressure.*/
        inline real S0(const real T);

        /*!\brief Entropy [J/(kg K)] at standard pressure.*/
        inline virtual real Ds0DT(const real T) const = 0;

        /*!\brief Molar entropy [J/(kmol K)] at standard pressure.*/
        inline real DS0DT(const real T);

        /*!\brief Entropy [J/(kg K)].*/
        inline virtual real s_p(const real pi, const real T) const;

        /*!\brief Entropy [J/(kg K)].*/
        inline virtual real s_rho(const real rhoi, const real T) const;

        /*!\brief The Gibbs free energy of species at standard pressure.*/
        inline virtual real g0(const real T) const;

        /*!\brief The standard-state molar Gibbs free energy.*/
        inline real G0(const real T) const;

        /*!\brief The Gibbs free energy of species at standard pressure.*/
        inline real Dg0DT(const real T) const;

        /*!\brief The standard-state molar Gibbs free energy.*/
        inline real DG0DT(const real T) const;

        /*!\brief The standard-state Helmholtz free energy.*/
        inline real a0(const real T) const;

        /*!\brief The standard-state molar Helmholtz free energy.*/
        inline real A0(const real T) const;

        /*!\brief Tha gas constant.*/
        inline real Ri() const;

        /*!\brief The molecular weight.*/
        inline real Wi() const;

        /**\brief Return the phase of the species.*/
        inline phaseType phase() const noexcept { return phase_; }

        /**\brief Return the phase of the species.*/
        inline void setPhase(const phaseType pT) { phase_ = pT; }
    };
} // namespace OpenHurricane

#include "thermo.inl"