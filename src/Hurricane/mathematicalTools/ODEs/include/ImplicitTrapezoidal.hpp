/*!
 * \file ImplicitTrapezoidal.hpp
 * \brief Headers of Implicit Trapezoidal Rule ODE solver.
 *        The subroutines and functions are in the <i>ImplicitTrapezoidal.cpp</i> file.
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
     * \brief The class of adaptive step Implicit Trapezoidal Rule ODE solver.
     */
    class ImplicitTrapezoidal : public ODEsSolver {
    private:
        realArray dfdt_;
        realArray dydt_;
        realSquareMatrix dfdy_;
        integerList pivotIndices_;

        realArray yk_;
        realArray yTemp_;
        realArray y05k_;

        hur_nodiscard real maxIncrement(const realArray &y0, const realArray &y,
                                        const realArray &e) const;

        bool isAutoUpdateJacobian_;

        integer stepToUpdateJac_;

        real ATolIncrement_;
        real RTolIncrement_;

        realArray dydtLast_;
        integer countStep_;

        realArray lastTimeStep_;

    public:
        declareClassNames;

        ImplicitTrapezoidal(const integer nEqs);

        ImplicitTrapezoidal(const integer nEqs, const real atol, const real rtol,
                            const integer maxStep = 1000);

        ImplicitTrapezoidal(const integer nEqs, const controller &cont);

        /**\brief Destructor.*/
        virtual ~ImplicitTrapezoidal() noexcept {}

        /*!\brief Return true if should reset the number of ODEs and reset it.*/
        virtual bool reset(const integer nEqsNew);

    public:
        /**\brief return maximum error.*/
        virtual real solve(const real t0, const realArray &y0, const realArray &dydt0,
                           const real dt, realArray &y);

        /*!\brief Only solve the ODEs without computing error.
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

        real NewtonIteration(const real t0, const realArray &y0, const realArray &dydt0,
                             const real dt, realArray &ykp1, bool &isConvergence);

        /*!\brief Solve the ODEs by giving the initial timestep
         *        and try to adjust the step according to the specified tolerance.
         * \param t - The time value.
         * \param dt - The initial timestep.
         * \param y - The solutions.
         */
        virtual void adaptiveSolve(real &t, real &dt, realArray &y);

    public:
        /*!\brief Solve the ODEs by giving the initial timestep in a period time
         *        and try to adjust the step according to the specified tolerance.
         * \param[in] t0 - The initial time value.
         * \param t0 - The end time value.
         * \param dt0 - The initial timestep.
         * \param y - The solutions.
         */
        virtual void solve(const real t0, real &tn, real &dt0, realArray &y);
    };

} // namespace OpenHurricane
