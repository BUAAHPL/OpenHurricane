/*!
 * \file eddyViscosity.hpp
 * \brief Header of eddy-viscosity hypothesis flow model.
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

#include "flowModel.hpp"
#include "viscousThermo.hpp"

namespace OpenHurricane {

    /*!\brief The class of eddy-viscosity hypothesis flow.*/
    class eddyViscosity : public flowModel {
    private:

        /*!\brief Schmidt number.*/
        real Sct_;

        /*!\brief Turbulent Prandtl number.*/
        real Prt_;

        /*!\brief The eddy viscosity field.*/
        cellRealArray mut_;

        /*!\brief Lower limit of mut*/
        real mutLow_;

        /*!\brief Higher limit of mut*/
        real mutHigh_;

    public:
        declareClassName(eddyViscosity);

        eddyViscosity(const runtimeMesh &mesh);

        eddyViscosity(const runtimeMesh &mesh, const controller &cont);

        /*!\brief Destructor.*/
        virtual ~eddyViscosity() noexcept {}

        /*!\brief Schmidt number.*/
        virtual hur_nodiscard real Sct() const noexcept { return Sct_; }

        /*!\brief Turbulent Prandtl number*/
        virtual hur_nodiscard real Prt() const noexcept { return Prt_; }
        // Member functions

        /*!\brief Return laminar dynamic viscosity field.*/
        virtual hur_nodiscard cellRealArray &mul() noexcept;

        /*!\brief Return laminar dynamic viscosity field.*/
        virtual hur_nodiscard const cellRealArray &mul() const noexcept;

        /*!\brief Return turbulent viscosity field.*/
        hur_nodiscard inline virtual cellRealArray &mut() noexcept { return mut_; }

        /*!\brief Return turbulent dynamic viscosity field.*/
        virtual hur_nodiscard const cellRealArray &mut() const noexcept;

        /*!\brief Return the effective dynamic viscosity field.
         *  The effective dynamic viscosity coefficient is defined as:
         * \f[\mu  = {\mu _l} + {\mu _t}\f]
         */
        virtual hur_nodiscard cellRealArray muEff();

        /*!\brief Return laminar thermal conductivity field.*/
        virtual hur_nodiscard cellRealArray &kappal() noexcept;

        /*!\brief Return laminar thermo conductivity field.*/
        virtual hur_nodiscard const cellRealArray &kappal() const noexcept;

        /*!\brief Return laminar thermal conductivity field.
         *  The turbulent thermal conductivity coefficient \f$\kappa_t\f$ is defined as:
         * \f[{\kappa _t} = {c_p}\frac{{{\mu _t}}}{{P{r_t}}}\f]
         */
        virtual hur_nodiscard cellRealArray &kappat() noexcept;

        /*!\brief Return the effective thermo conductivity field.
         *  The effective thermal conductivity is defined as:
         * \f[\kappa  = {\kappa _l} + {\kappa _t}\f]
         */
        virtual hur_nodiscard cellRealArray kappaEff();

        /*!\brief Return the effective mu/Pr:.
         *  \f[\frac{{{\mu _l}}}{{P{r_l}}} + \frac{{{\mu _t}}}{{P{r_t}}}\f]
         */
        virtual hur_nodiscard cellRealArray keEff() const;

        /*!\brief Return the effective mu/Pr of cell <i>i</i>.
         *  \f[\frac{{{\mu _l}}}{{P{r_l}}} + \frac{{{\mu _t}}}{{P{r_t}}}\f]
         */
        virtual hur_nodiscard real keEff(const integer cellI) const;

        /*!\brief Lower limit of mut*/
        inline virtual hur_nodiscard real mutLow() const noexcept { return mutLow_; }

        /*!\brief Higher limit of mut*/
        inline virtual hur_nodiscard real mutHigh() const noexcept { return mutHigh_; }

        inline virtual hur_nodiscard string typeName() const { return className_; }
    };
} // namespace OpenHurricane