/*!
 * \file AUSMPlusUP.hpp
 * \brief Headers of AUSM+ -up scheme.
 *        The subroutines and functions are in the <i>AUSMPlusUP.cpp</i> file.
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
#include "upwindScheme.hpp"

namespace OpenHurricane {

    /*!
     * \brief The class of AUSMPlusUP scheme.
     */
    class AUSMPlusUP : public upwindScheme {
    private:
        real Kp_;
        real Ku_;
        real sigma_;
        real beta_;

        hur_nodiscard inline real alpha(const real fa) const;

        hur_nodiscard inline real fa(const real Mo) const;

        hur_nodiscard inline real M1P(const real Ma) const;
        hur_nodiscard inline real M1M(const real Ma) const;

        hur_nodiscard inline real M2P(const real Ma) const;
        hur_nodiscard inline real M2M(const real Ma) const;

        hur_nodiscard inline real M4P(const real Ma) const;
        hur_nodiscard inline real M4M(const real Ma) const;

        hur_nodiscard inline real P5P(const real Ma, const real _alpha) const;
        hur_nodiscard inline real P5M(const real Ma, const real _alpha) const;

        mutable realArray fp_;

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
        inline AUSMPlusUP() = delete;

        AUSMPlusUP(const controller &cont, const runtimeMesh &mesh, flowModel &flow);

        /**\brief Destructor.*/
        virtual ~AUSMPlusUP() noexcept {}
    };

} // namespace OpenHurricane

#include "AUSMPlusUP.inl"