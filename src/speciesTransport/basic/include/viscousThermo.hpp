/*!
 * \file viscousThermo.hpp
 * \brief Headers of class of viscous thermo.
 *        The subroutines and functions are in the <i>viscousThermo.cpp</i> file.
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
#include "rhoThermo.hpp"

namespace OpenHurricane {
    class viscousThermo : public rhoThermo {
    private:
        /*!\brief Dynamic viscosity field.*/
        mutable cellRealArray mu_;

        /*!\brief Laminar thermo conductivity*/
        mutable cellRealArray kappa_;

    public:

        /*!\brief Construct from mesh*/
        viscousThermo(const runtimeMesh &mesh);

        /*!\brief Construct from mesh*/
        viscousThermo(const runtimeMesh &mesh, const controller &cont);

        /*!\brief Destructor.*/
        virtual ~viscousThermo() noexcept {}

        /*!\brief Return const access to dynamic viscosity field.*/
        virtual hur_nodiscard const cellRealArray &mu() const noexcept;

        /*!\brief Return access to dynamic viscosity field.*/
        virtual hur_nodiscard cellRealArray &mu() noexcept;

        /*!\brief Return const access to Laminar thermo conductivity.*/
        virtual hur_nodiscard const cellRealArray &kappa() const noexcept;

        /*!\brief Return access to Laminar thermo conductivity.*/
        virtual hur_nodiscard cellRealArray &kappa() noexcept;
    };
} // namespace OpenHurricane