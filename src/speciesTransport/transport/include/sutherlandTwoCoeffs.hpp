﻿/*!
 * \file sutherlandTwoCoeffs.hpp
 * \brief Header of transport properties by sutherland two coeffs.
 *       The subroutines and functions are in the <i>sutherlandTwoCoeffs.inl</i> file.
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
#include "transport.hpp"

namespace OpenHurricane {
    /*!\brief The class of sutherland two coeffs transport properties.*/
    class sutherlandTwoCoeffs : public transport {
    private:
        real C1_;

        real C2_;

    public:
        declareClassNames;

        inline sutherlandTwoCoeffs(const speciesList &sp, const integer index, const real C1,
                                   const real C2, const real Prl);

        inline sutherlandTwoCoeffs(const speciesList &sp, const integer index,
                                   const controller &cont);

        sutherlandTwoCoeffs(const sutherlandTwoCoeffs &tra);

        sutherlandTwoCoeffs(const sutherlandTwoCoeffs &tra, const speciesList &sp);

        /*!\brief Return a clone.*/
        virtual hur_nodiscard uniquePtr<transport> clone() const {
            return uniquePtr<transport>(new sutherlandTwoCoeffs(*this));
        }

        /*!\brief Return a clone by given new species list.*/
        virtual hur_nodiscard uniquePtr<transport> clone(const speciesList &sp) const {
            return uniquePtr<transport>(new sutherlandTwoCoeffs(*this, sp));
        }

        /*!\brief Destructor.*/
        virtual ~sutherlandTwoCoeffs() noexcept;

        /*!
         * \brief Return the dynamic viscosity [kg/(m s)].
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         */
        hur_nodiscard inline virtual real mu(const real p, const real T) const;

        /*!
         * \brief Return the thermal conductivity [W/mK].
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard inline virtual real kappa(const real p, const real T, const real cpi) const;

        /*!
         * \brief Return the thermal conductivity [W/mK].
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] mui - Dynamic viscosity [kg/(m s)].
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard inline virtual real kappa(const real p, const real T, const real mui,
                                                const real cpi) const;

        sutherlandTwoCoeffs &operator=(const sutherlandTwoCoeffs &);
    };
} // namespace OpenHurricane

#include "sutherlandTwoCoeffs.inl"