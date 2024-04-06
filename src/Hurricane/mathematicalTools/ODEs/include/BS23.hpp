/*!
 * \file BS23.hpp
 * \brief Headers of Bogacki-Shampine 2/3 ODEs embedded pair solver.
 *        The subroutines and functions are in the <i>BS23.cpp</i> file.
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
     * \brief The class of Bogacki-Shampine 2/3 ODEs embedded pair solver.
     */
    class BS23 : public ODEsSolver {
    private:
        realArray s2_;
        realArray s3_;
        realArray s4_;
        realArray error_;

    public:
        declareClassNames;

        BS23(const integer nEqs);

        BS23(const integer nEqs, const real atol, const real rtol, const integer maxStep = 1000);

        BS23(const integer nEqs, const controller &cont);

        virtual ~BS23() noexcept {}

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
    };

} // namespace OpenHurricane
