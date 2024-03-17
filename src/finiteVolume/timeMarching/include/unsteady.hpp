/*!
 * \file unsteady.hpp
 * \brief Headers of class of the unsteady time schemes.
 *        The subroutines and functions are in the <i>unsteady.cpp</i> file.
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
#include "timeMarching.hpp"

namespace OpenHurricane {

    /*!\brief The base class for unsteady time scheme.*/
    class unsteady {
    private:
        /*!\brief Hold reference to the runtime mesh.*/
        const runtimeMesh &mesh_;

        /*!\brief Thermo.*/
        const flowModel &flowM_;

        solver &solverM_;

    public:
        declareClassName(unsteady);

        /*!\brief Disallow null constructor.*/
        unsteady() = delete;

        /*!\brief Disallow copy constructor.*/
        unsteady(const unsteady &) = delete;

        unsteady(const runtimeMesh &mesh, const flowModel &flowM, solver &_solver);

        /*!\brief Destructor.*/
        virtual ~unsteady() noexcept {}

        void getTimeStep(realArray &dt) const;

        void getShockFactor(cellRealArray &shockFactor_) const;
    };

} // namespace OpenHurricane