/*!
 * \file constCp.hpp
 * \brief Header of thermo properties by constCp.
 *       The subroutines and functions are in the <i>constCp.cpp</i> file.
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

    class constCp : public thermo {
    private:
        /*!\brief Heat capacity at constant pressure [J/(kg K)].*/
        real cp_;

        /*!\brief Chemical enthalpy [J/kg].*/
        real hc_;

    public:
        declareClassNames;

        constCp(const controller &cont, const equationOfState &st, const integer id);

        constCp(const equationOfState &st, const real cp, const real hc, const integer id);

        inline constCp(const constCp &);
        constCp &operator=(const constCp &cC);

        inline constCp(const constCp &, const string &);

        virtual hur_nodiscard uniquePtr<thermo> clone() const {
            return uniquePtr<thermo>(new constCp(*this));
        }

        virtual ~constCp() noexcept {}

        /*!\brief Heat capacity at constant pressure [J/(kg K)] at standard pressure (1 atm).*/
        hur_nodiscard inline virtual real cp0(const real T) const noexcept;

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

#include "constCp.inl"