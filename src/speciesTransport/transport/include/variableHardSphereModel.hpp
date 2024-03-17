/*!
 * \file variableHardSphereModel.hpp
 * \brief Header of transport properties by the variable hard sphere model.
 *       The subroutines and functions are in the <i>variableHardSphereModel.inl</i> file.
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
    /*!\brief The class of the variable hard sphere model.*/
    class variableHardSphereModel : public transport {
    private:
        /** \brief VHS temperature exponent. */
        real omega_;

        /** \brief Reference temperature [K]. */
        real Tref_;

        /** \brief Reference molecular diameter [m]. */
        real dref_;

        /** \brief The mass of one molecule. */
        real m_;

        /** \brief The dynamic viscosity [kg/(m s)] of reference state. */
        real muref_;

        void calcMuRef();

    public:
        declareClassNames;

        inline variableHardSphereModel(const speciesList &sp, const integer index, const real omega,
                                       const real Tref, const real dref, const real Prl);

        inline variableHardSphereModel(const speciesList &sp, const integer index,
                                       const controller &cont);

        variableHardSphereModel(const variableHardSphereModel &tra);

        variableHardSphereModel(const variableHardSphereModel &tra, const speciesList &sp);

        /*!\brief Return a clone.*/
        virtual hur_nodiscard uniquePtr<transport> clone() const {
            return uniquePtr<transport>(new variableHardSphereModel(*this));
        }

        /*!\brief Return a clone by given new species list.*/
        virtual hur_nodiscard uniquePtr<transport> clone(const speciesList &sp) const {
            return uniquePtr<transport>(new variableHardSphereModel(*this, sp));
        }

        /*!\brief Destructor.*/
        virtual ~variableHardSphereModel() noexcept;

        /*!
         * \brief Return the dynamic viscosity [kg/(m s)].
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         */
        inline virtual real mu(const real p, const real T) const;

        /*!
         * \brief Return the thermal conductivity [W/mK].
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        inline virtual real kappa(const real p, const real T, const real cpi) const;

        /*!
         * \brief Return the thermal conductivity [W/mK].
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] mui - Dynamic viscosity [kg/(m s)].
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        inline virtual real kappa(const real p, const real T, const real mui, const real cpi) const;

        variableHardSphereModel &operator=(const variableHardSphereModel &);
    };
} // namespace OpenHurricane

#include "variableHardSphereModel.inl"