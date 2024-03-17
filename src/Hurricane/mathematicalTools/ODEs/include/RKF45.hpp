/*!
 * \file RKF45.hpp
 * \brief Headers of 4/5th Order Runge-Kutta-Fehlberg embedded pair ODE solver.
 *        The subroutines and functions are in the <i>RKF45.cpp</i> file.
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

#include "ODEsSolver.hpp"

namespace OpenHurricane {
    /**
     * \brief The class of 4/5th Order Runge-Kutta-Fehlberg embedded pair ODE solver.
     */
    class RKF45 : public ODEsSolver {
    private:

        static constexpr real c2 = 1.0 / 4.0;
        static constexpr real c3 = 3.0 / 8.0;
        static constexpr real c4 = 12.0 / 13.0;
        static constexpr real c5 = 1.0;
        static constexpr real c6 = 1.0 / 2.0;

        static constexpr real a21 = 1.0 / 4.0;
        static constexpr real a31 = 3.0 / 32.0;
        static constexpr real a32 = 9.0 / 32.0;
        static constexpr real a41 = 1932.0 / 2197.0;
        static constexpr real a42 = -7200.0 / 2197.0;
        static constexpr real a43 = 7296.0 / 2197.0;
        static constexpr real a51 = 439.0 / 216.0;
        static constexpr real a52 = -8.0;
        static constexpr real a53 = 3680.0 / 513.0;
        static constexpr real a54 = -845.0 / 4104.0;
        static constexpr real a61 = -8.0 / 27.0;
        static constexpr real a62 = 2.0;
        static constexpr real a63 = -3544.0 / 2565.0;
        static constexpr real a64 = 1859.0 / 4104.0;
        static constexpr real a65 = -11.0 / 40.0;

        static constexpr real b1 = 16.0 / 135.0;
        static constexpr real b3 = 6656.0 / 12825.0;
        static constexpr real b4 = 28561.0 / 56430.0;
        static constexpr real b5 = -9.0 / 50.0;
        static constexpr real b6 = 2.0 / 55.0;

        static constexpr real e1 = 25.0 / 216.0 - b1;
        static constexpr real e3 = 1408.0 / 2565.0 - b3;
        static constexpr real e4 = 2197.0 / 4104.0 - b4;
        static constexpr real e5 = -1.0 / 5.0 - b5;
        static constexpr real e6 = -b6;

        /**\brief Temporary fields.*/
        realArray s2_;
        realArray s3_;
        realArray s4_;
        realArray s5_;
        realArray s6_;
        realArray error_;
        realArray yTmp_;

    public:
        declareClassNames;

        RKF45(const integer nEqs);

        RKF45(const integer nEqs, const real atol, const real rtol, const integer maxStep = 1000);

        RKF45(const integer nEqs, const controller &cont);

        virtual ~RKF45() noexcept {}

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
