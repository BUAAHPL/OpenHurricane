/*!
 * \file Euler.hpp
 * \brief Headers of Euler ODEs solver.
 *        The subroutines and functions are in the <i>Euler.cpp</i> file.
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

#include "ODEsSolver.hpp"

namespace OpenHurricane {

    /**
     * \brief The class of adaptive step Euler solver.
     */
    class Euler : public ODEsSolver {

    private:
        realArray yTemph2_;
        realArray error_;
        realArray dydt2_;

    public:
        declareClassNames;

        Euler(const integer nEqs);

        Euler(const integer nEqs, const real atol, const real rtol, const integer maxStep = 1000);

        Euler(const integer nEqs, const controller &cont);

        virtual ~Euler() noexcept {}

        /*!\brief Return true if should reset the number of ODEs and reset it.*/
        virtual bool reset(const integer nEqsNew);

    public:
        /**\brief return maximum error.*/
        virtual real solve(const real t0, const realArray &y0, const realArray &dydt0,
                           const real dt, realArray &y);

        /*!\brief Only solve the ODEs without computing error only for explicit solver.
         * \param[in] t0 - The initial time.
         * \param[in] y0 - The initial states.
         * \param[in] dydt0 - The derivatives.
         * \param[in] dt - The timestep.
         * \param[out] y - The solutions.
         */
        virtual void solveOnly(const real t0, const realArray &y0, const realArray &dydt0,
                               const real dt, realArray &y);

    protected:
        virtual real solve(const real t0, const realArray &y0, const realArray &dydt0,
                           const real dt, const integerListList &Gml, const integer im,
                           realArray &y);

        /*!\brief Get error estimate.
         * \param[in] t0 - The initial time.
         * \param[in] tbase - The base physical time, and the time the solver should reach.
         * \param[out] tReached - The time reached by solver.
         * \param[in] dt - The character time scale array.
         * \param[in] Gml - The groups of different time scale.
         * \param y - The initial values and the solutions.
         */
        virtual real getError(const real t0, const real tbase, real &tReached, const realArray &dt,
                              const integerListList &Gml, realArray &y);
    };

} // namespace OpenHurricane
