/*!
 * \file expCFL.hpp
 * \brief Header of exponential law adapt CFL number.
 *       The subroutines and functions are in the <i>expCFL.cpp</i> file.
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

#include "CFL.hpp"
#include "cellArrays.hpp"
#include "iteration.hpp"

namespace OpenHurricane {
    /*!\brief The class of exponential law adapt CFL number.*/
    class expCFL : public CFL {
    private:
        integer CFLConst_;

        real getCFLFromStep(const integer cstep, const bool isPrintCFL = false) const;

    public:
        declareClassNames;

        /*!\brief Disallow null constructor.*/
        expCFL() = delete;

        /*!\brief Construct from controller.*/
        expCFL(const iteration &iter, const runtimeMesh &mesh, const controller &cont);

        /*!\brief Disallow copy constructor.*/
        expCFL(const expCFL &cfl) = delete;

        /*!\brief Destructor.*/
        virtual ~expCFL() noexcept {}

        virtual integer CFLConst() const noexcept;

        /*!\brief Return the CFL number for current step.
         * \param[in] cstep - The current step.
         * \return The CFL number for current step.
         */
        virtual real getCFL() const;
    };
} // namespace OpenHurricane
