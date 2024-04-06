/*!
 * \file laminarFlow.hpp
 * \brief Header of laminar flow model.
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

#include "flowModel.hpp"

#include "viscousThermo.hpp"

namespace OpenHurricane {
    /*!\brief The class of laminar flow.*/
    class laminarFlow : public flowModel {
    public:
        declareClassName(laminarFlow);

        laminarFlow(const runtimeMesh &mesh);

        laminarFlow(const runtimeMesh &mesh, const controller &cont);

        /*!\brief Destructor.*/
        virtual ~laminarFlow() noexcept {}


        /*!\brief Return laminar dynamic viscosity field.*/
        virtual hur_nodiscard cellRealArray &mul() noexcept;

        /*!\brief Return laminar dynamic viscosity field.*/
        virtual hur_nodiscard const cellRealArray &mul() const noexcept;

        /*!\briefReturn the effective dynamic viscosity field.*/
        virtual hur_nodiscard cellRealArray muEff();

        /*!\brief Return laminar thermo conductivity field.*/
        virtual hur_nodiscard cellRealArray &kappal() noexcept;

        /*!\brief Return laminar thermo conductivity field.*/
        virtual hur_nodiscard const cellRealArray &kappal() const noexcept;

        /*!\briefReturn the effective thermo conductivity field.*/
        virtual hur_nodiscard cellRealArray kappaEff();

        /*!\briefReturn the effective mu/Pr.*/
        virtual hur_nodiscard cellRealArray keEff() const;

        /*!\briefReturn the effective mu/Pr.*/
        virtual hur_nodiscard real keEff(const integer cellI) const;

        inline virtual hur_nodiscard string typeName() const { return className_; }
    };
} // namespace OpenHurricane