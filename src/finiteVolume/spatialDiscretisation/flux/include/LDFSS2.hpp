/*!
 * \file LDFSS2.hpp
 * \brief Headers of Low-Diffusion Flux Splitting Scheme.
 *        The subroutines and functions are in the <i>LDFSS2.cpp</i> file.
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
#include "upwindScheme.hpp"

namespace OpenHurricane {

    /*!
     * \brief The class of LDFSS2 scheme.
     */
    class LDFSS2 : public upwindScheme {
    private:
        mutable realArray pFlux_;
        mutable realArray massFluxL_;
        mutable realArray massFluxR_;

    public:
        static const real delta_;

    public:
        /**
         * \brief The function for calculating convective flux.
         * \param[in] rhol - The left state density.
         * \param[in] rhor - The right state density.
         * \param[in] VL - The left state velocity.
         * \param[in] VR - The right state velocity.
         * \param[in] pl - The left state pressure.
         * \param[in] pr - The right state pressure.
         * \param[in] gl - The left state specific heat ratio.
         * \param[in] gr - The right state specific heat ratio.
         * \param[in] el - The left state total energy.
         * \param[in] er - The right state total energy.
         * \param[in] cl - The left state acoustic speed.
         * \param[in] cr - The right state acoustic speed.
         * \param[in] faceArea - The area vector of this face.
         * \param[in] MasInf - The Mach number of infinity.
         * \param[out] flux - The flux of this face.
         * \return The reaction rate conatants.
         */
        virtual void calcFlux(const real rhol, const real rhor, const vector &VL, const vector &VR,
                              const real pl, const real pr, const real gl, const real gr,
                              const real el, const real er, const real cl, const real cr,
                              const vector &faceArea, const real blend, realArray &flux) const;

    public:
        declareClassNames;

        /**\brief Construct null.*/
        LDFSS2() = delete;

        LDFSS2(const controller &cont, const runtimeMesh &mesh, flowModel &flow);

        /**\brief Destructor.*/
        virtual ~LDFSS2() noexcept {}

    private:
        hur_nodiscard inline real M1P(const real Ma) const { return 0.5 * (Ma + mag(Ma)); }

        hur_nodiscard inline real M1M(const real Ma) const { return 0.5 * (Ma - mag(Ma)); }

        hur_nodiscard inline real M2P(const real Ma) const { return 0.25 * sqr(Ma + real(1.0)); }

        hur_nodiscard inline real M2M(const real Ma) const { return -0.25 * sqr(Ma - real(1.0)); }

        hur_nodiscard inline real P5P(const real Ma, const real _alpha) const {
            real p;
            if (std::abs(Ma) > 1.0) {
                p = M1P(Ma) / Ma;
            } else {
                p = M2P(Ma) * ((2.0 - Ma) - 16.0 * _alpha * Ma * M2M(Ma));
            }
            return p;
        }

        hur_nodiscard inline real P5M(const real Ma, const real _alpha) const {
            real p;
            if (std::abs(Ma) > 1.0) {
                p = M1M(Ma) / Ma;
            } else {
                p = M2M(Ma) * ((-2.0 - Ma) + 16.0 * _alpha * Ma * M2P(Ma));
            }
            return p;
        }
    };

} // namespace OpenHurricane