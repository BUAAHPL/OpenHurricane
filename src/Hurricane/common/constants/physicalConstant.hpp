/*!
 * \file physicalConstant.hpp
 * \brief All the information about the definition of the physical constants.
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

#include "mathConstants.hpp"
#include "real.hpp"

namespace OpenHurricane {
    namespace constant {
        /**
         * \brief The namespace of the physical constants.
         */
        namespace physicalConstant {
            /**
             *\brief Standard acceleration of gravity (gee, free-fall on Earth)
             *       (Units: [m/s2]).
             */
            constexpr hur_inline_var real gn = 9.80665;

            /**\brief Speed of light in vacuum (Unit: [m/s]).*/
            constexpr hur_inline_var real c = 299792458;

            /**\brief Planck constant (Unit: [J.s]).*/
            constexpr hur_inline_var real h = 6.62607015e-34;

            /**\brief Reduced Planck constant (Unit: [J/s]).*/
            constexpr hur_inline_var real hr = h / mathConstants::twoPi;

            /**\brief Newtonian constant of gravitation (Unit: [N.m2/kg2]).*/
            constexpr hur_inline_var real G = 6.67408e-11;

            /**\brief Elementary charge (Unit: [C]).*/
            constexpr hur_inline_var real e = 1.6021766208e-19;

            /**\brief Magnetic permeability of free space (Unit: [H/m]).*/
            constexpr hur_inline_var real mu0 = 4.0 * mathConstants::pi * 1e-7;

            /**\brief Electric constant (vacuum permittivity) (Unit: [F/m]).*/
            constexpr hur_inline_var real epsilon0 = 1.0 / (mu0 * sqr(c));

            /**\brief Characteristic impedance of a vacuum (Unit: [ohm]).*/
            constexpr hur_inline_var real Z0 = mu0 * c;

            /**\brief Coulomb constant(Units: [N.m2/C2]).*/
            constexpr hur_inline_var real kappa = 1.0 / (4.0 * mathConstants::pi * epsilon0);

            /**\brief Conductance quantum (Units: [S]).*/
            constexpr hur_inline_var real G0 = (2.0 * sqr(e)) / h;

            /**\brief Josephson constant (Units: [Hz/V]).*/
            constexpr hur_inline_var real KJ = 2.0 * e / h;

            /**\brief Magnetic flux quantum (Units: [Wb]).*/
            constexpr hur_inline_var real phi0 = h / (2.0 * e);

            /**\brief Von Klitzing constant (Units: [ohm]).*/
            constexpr hur_inline_var real RK = h / sqr(e);

            /**\brief Electron mass (Unit: [kg]).*/
            constexpr hur_inline_var real me = 9.10938356e-31;

            /**\brief Bohr magneton (Units: []).*/
            constexpr hur_inline_var real muB = e * hr / (2.0 * me);

            /**\brief Fine-structure constant (Units: []).*/
            constexpr hur_inline_var real alpha = sqr(e) / (2.0 * epsilon0 * h * c);

            /**\brief Rydberg constant (Units: [1/m]).*/
            constexpr hur_inline_var real Rinf = sqr(alpha) * me * c / (2.0 * h);

            /**\brief Bohr radius (Units: [m]).*/
            constexpr hur_inline_var real a0 = alpha / (4.0 * mathConstants::pi * Rinf);

            /**\brief Classical electron radius (Units: [m]).*/
            constexpr hur_inline_var real re =
                sqr(e) / (4.0 * mathConstants::pi * epsilon0 * me * sqr(c));

            /**\brief Hartree energy (Units: [J]).*/
            constexpr hur_inline_var real Eh = 2.0 * Rinf * h * c;

            /**\brief Proton mass (Units: [kg]).*/
            constexpr hur_inline_var real mp = 1.672621898e-27;

            /**\brief Electronvolt (Units: [J])
             *		1eV = e * 1V
             *		1V = 1 J/C
             *		e is elementary charge.*/
            constexpr hur_inline_var real eV = e;

            /**\brief Thermochemical calorie (Units: [J]).*/
            constexpr hur_inline_var real calTh = 4.184;

            /**\brief 4 ¡ãC calorie (Units: [J]).*/
            constexpr hur_inline_var real cal4 = 4.204;

            /**\brief 15 ¡ãC calorie (Units: [J]).*/
            constexpr hur_inline_var real cal15 = 4.1855;

            /**\brief 20 ¡ãC calorie (Units: [J]).*/
            constexpr hur_inline_var real cal20 = 4.182;

            /**\brief Mean calorie (Units: [J]).*/
            constexpr hur_inline_var real calMean = 4.190;

            /**\brief International Steam table calorie (1929) (Units: [J]).*/
            constexpr hur_inline_var real calITold = 4.1868;

            /**\brief International Steam table calorie (1956) (Units: [J]).*/
            constexpr hur_inline_var real calIT = 4.1868;

            /**\brief Default calorie unit (Units: [J])
             *		equal to thermochemical calorie:
             *		4.184 J.*/
            constexpr hur_inline_var real cal = calTh;

            /**\brief Kilowatt-hour (Unit: [J]).*/
            constexpr hur_inline_var real kWh = 3600.0 * 1e3;

            /**\brief Avogadro constant (Unit: [mol-1]).*/
            constexpr hur_inline_var real NA = 6.02214076e23;

            /**\brief Boltzmann constant (Unit: [J/K]).*/
            constexpr hur_inline_var real k = 1.380649e-23;

            /**\brief Atomic mass unit (Unit: [kg]).*/
            constexpr hur_inline_var real AMU = 1.6605402e-27;

            /**\brief Faraday constant (Unit: [C/mol]).*/
            constexpr hur_inline_var real F = e * NA;

            /**\brief Universal gas constant (Unit: [J/mol/K]).*/
            constexpr hur_inline_var real R = k * NA;

            /**\brief Standard molar volume of gas (Unit: [L/mol])
             *	at T=273.15K and p=101.325kPa.
             */
            constexpr hur_inline_var real Vm = 22.414110;

            /**\brief Stefan-Boltzmann constant (Unit: [W/m2/K4]).*/
            constexpr hur_inline_var real sigma =
                (sqr(mathConstants::pi) / 60.0) * pow4(k) / (OpenHurricane::pow3(hr) * sqr(c));

            /**\brief Wien displacement law constant(Unit: [m.K]).*/
            constexpr hur_inline_var real b = (h * c / k) / 4.965114231;
            //4.965114231 is the solution of the equation: e^-x +x/5 = 1
            //const real b = 0.002897;

            /**\brief First radiation constant(Unit: [W/m2]).*/
            constexpr hur_inline_var real c1 = mathConstants::twoPi * h * sqr(c);
            //const real c1 = 3.7417749e-16;

            /**\brief Second radiation constant(Unit: [m.K]).*/
            constexpr hur_inline_var real c2 = h * c / k;
            //const real c2 = 0.01438769;

            /**\brief Universal gas constant (Unit: [J/(kmol K)]).*/
            constexpr hur_inline_var real Ru = R * 1e3;

            /**\brief Standard atmosphere pressure (Unit: [Pa]).*/
            constexpr hur_inline_var real Patm = 1.01325e5;

            /**\brief Standard temperature (Unit: [K]).*/
            constexpr hur_inline_var real Tstd = 273.15;
        } // namespace physicalConstant
    }     // namespace constant
} // namespace OpenHurricane