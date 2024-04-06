/*!
 * \file transport.hpp
 * \brief Header of transport properties.
 *       The subroutines and functions are in the <i>transport.cpp</i> file.
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

#include "commonInclude.hpp"
#include "logFile.hpp"
#include "objectFactory.hpp"
#include "realArray.hpp"
#include "speciesList.hpp"

#include "smartPointerList.hpp"

namespace OpenHurricane {
    /*!\brief The base class of transport properties.*/
    class transport {
    private:
        /*!\brief Const reference to species list.*/
        const speciesList &species_;

        /*!\brief Index of species.*/
        integer index_;

        /*!\brief Laminar Prandtl number.*/
        real Pr_;

    public:
        declareClassNames;

        declareObjFty(transport, controller,
                      (const speciesList &sp, const integer index, const controller &cont),
                      (sp, index, cont));

        transport(const speciesList &sp, const integer index, const real Prl);

        transport(const speciesList &sp, const integer index, const controller &cont);

        transport(const transport &tra);

        /*!\brief Construct as copy and given new species list.*/
        transport(const transport &tra, const speciesList &sp);

        /*!\brief Return a clone.*/
        virtual hur_nodiscard uniquePtr<transport> clone() const = 0;

        /*!\brief Return a clone by given new species list.*/
        virtual hur_nodiscard uniquePtr<transport> clone(const speciesList &sp) const = 0;

        /*!\brief Select null constructed.*/
        static uniquePtr<transport> creator(const speciesList &sp, const integer index,
                                            const controller &cont);

        /*!\brief Destructor.*/
        virtual ~transport() noexcept {}

        /*!\brief Return const access to the species list.*/
        hur_nodiscard inline const speciesList &species() const noexcept;

        hur_nodiscard inline integer index() const noexcept;

        /*!\brief Return laminar Prandtl number.*/
        hur_nodiscard inline real Pr() const noexcept;

        /*!\brief Return laminar Prandtl number.*/
        inline real setPr(const real newPr) noexcept{
            real oldPr = Pr_;
            Pr_ = newPr;
            return oldPr;
        }

        /*!\brief The Lennard-Jones collision diameter in angstroms.*/
        hur_nodiscard inline virtual real sigma() const noexcept { return 0.0; }

        // Transport properties

        /*!
         * \brief Return the dynamic viscosity [kg/(m s)].
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         */
        hur_nodiscard virtual real mu(const real p, const real T) const = 0;

        /*!
         * \brief Return the thermal conductivity [W/mK].
         * \param[in] p - Static pressure.
         * \param[in] T - Static temperature.
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard virtual real kappa(const real p, const real T, const real cpi) const = 0;

        /*!
         * \brief Return the thermal conductivity [W/mK].
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] mui - Dynamic viscosity [kg/(m s)].
         * \param[in] cpi - Mass-based heat capacity at constant pressure [J/(kg K)].
         */
        hur_nodiscard virtual real kappa(const real p, const real T, const real mui,
                                         const real cpi) const = 0;

        /*!
         * \brief Return the Species mass diffusivity.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] xi - Molar fraction of species.
         */
        hur_nodiscard inline virtual real Di(const real p, const real T, const realArray &xi,
                                             const PtrList<transport> &tranPtr) const;

        /*!
         * \brief Return the Species mass diffusivity.
         * \param[in] p - Static pressure [Pa].
         * \param[in] T - Static temperature [K].
         * \param[in] xi - Molar fraction of species.
         */
        hur_nodiscard inline virtual real Di(const real p, const real T, const realArray &xi,
                                             const PtrList<transport> &tranPtr,
                                             realArray &Dij) const;

        hur_nodiscard virtual real Di(const real p, const real T, const real Tv1p5,
                                      const realArray &xi, const PtrList<transport> &tranPtr,
                                      realArray &Dij) const;

        hur_nodiscard virtual real Di(const real p, const real T, const real Tv1p5,
                                      const realArray &xi, const PtrList<transport> &tranPtr) const;

        transport &operator=(const transport &t);
    };
} // namespace OpenHurricane

#include "transport.inl"