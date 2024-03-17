/*!
 * \file piecewiseLinear.hpp
 * \brief Header of thermo properties by piecewise-linear.
 *       The subroutines and functions are in the <i>piecewiseLinear.cpp</i> file.
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

#include "FixedList.hpp"
#include "commonInclude.hpp"
#include "controller.hpp"
#include "real.hpp"
#include "thermo.hpp"

namespace OpenHurricane {

    class piecewiseLinear : public thermo {
    private:
        /*!\brief Heat capacity at constant pressure [J/(kg K)] at different temperature.*/
        realArray cp_;

        /** \brief Temperature array. */
        realArray T_;

        /*!\brief Chemical enthalpy [J/kg].*/
        real hc_;

    public:
        declareClassNames;

        piecewiseLinear(const controller &cont, const equationOfState &st, const integer id);

        piecewiseLinear(const equationOfState &st, const realArray cpl, const realArray Tl,
                        const real hc, const integer id);

        inline piecewiseLinear(const piecewiseLinear &);
        piecewiseLinear &operator=(const piecewiseLinear &cC);

        inline piecewiseLinear(const piecewiseLinear &, const string &);

        /*!\brief Return a clone.*/
        virtual hur_nodiscard uniquePtr<thermo> clone() const {
            return uniquePtr<thermo>(new piecewiseLinear(*this));
        }

        virtual ~piecewiseLinear() noexcept {}

        /*!\brief Heat capacity at constant pressure [J/(kg K)] at standard pressure (1 atm).*/
        hur_nodiscard inline virtual real cp0(const real T) const noexcept;

        hur_nodiscard inline virtual real inteCp0dT(const real T1, const real T2) const noexcept;

        /*!\brief Absolute Enthalpy [J/kg] at standard pressure.*/
        hur_nodiscard inline virtual real ha0(const real T) const noexcept;

        /*!\brief Absolute Enthalpy [J/kg] at standard pressure.*/
        hur_nodiscard inline virtual real Dha0DT(const real T) const noexcept;

        /*!\brief Chemical enthalpy [J/kg].*/
        hur_nodiscard inline virtual real hc() const noexcept;

        /*!\brief Entropy [J/(kg K)] at standard pressure.*/
        hur_nodiscard inline virtual real s0(const real T) const;

        /*!\brief Entropy [J/(kg K)] at standard pressure.*/
        hur_nodiscard inline virtual real Ds0DT(const real T) const noexcept;
    };

} // namespace OpenHurricane

#include "piecewiseLinear.inl"