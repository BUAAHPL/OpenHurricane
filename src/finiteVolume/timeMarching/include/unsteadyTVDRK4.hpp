/*!
 * \file unsteadyTVDRK4.hpp
 * \brief Headers of class of the unsteady 4th order TVD Runge-Kutta.
 *        The subroutines and functions are in the <i>unsteadyTVDRK4.cpp</i> file.
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
#include "TVDRK4.hpp"
#include "unsteady.hpp"

namespace OpenHurricane {

    /*!\brief The base class for unsteady TVD RK scheme.*/
    class unsteadyTVDRK4 : public TVDRK4, public unsteady {
    public:
        declareClassName(unsteadyTVDRK4);

        /*!\brief Disallow null constructor.*/
        unsteadyTVDRK4() = delete;

        /*!\brief Disallow copy constructor.*/
        unsteadyTVDRK4(const unsteadyTVDRK4 &) = delete;

        unsteadyTVDRK4(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
                       solver &_solver, cellVectorArray &v);

        /*!\brief Destructor.*/
        virtual ~unsteadyTVDRK4() noexcept {}

        virtual void timeStep();
    };

} // namespace OpenHurricane