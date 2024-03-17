/*!
 * \file referenceValues.hpp
 * \brief Header of reference values
 *       The subroutines and functions are in the <i>referenceValues.cpp</i> file.
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

#include "dataStructure.hpp"
#include "realArray.hpp"

namespace OpenHurricane {
    class controller;

    class referenceValues {
    private:
        /** \brief Mach Number. */
        real Ma_;

        /** \brief Velocity Magnitude [m/s]. */
        real u_;

        /** \brief Static pressure [Pascal]*/
        real p_;

        /** \brief Static temperature. [K] */
        real T_;

        /** \brief Density [kg/m^3] */
        real rho_;

        /** \brief Viscosity [kg/(m s)] */
        real mu_;

        /** \brief Ratio of specific heats. */
        real gamma_;

        /** \brief The reference area, which is used to compute the force and the moment coefficients. [m^2] */
        real area_;

        /** \brief The reference length, which is used to compute the moment coefficient. [m] */
        real length_;

        /** \brief The reference stagnation temperature [K] */
        real Tt_;

        /** \brief Gas constant. */
        real R_;

    public:
        /** \brief Null constructor. */
        inline referenceValues()
            : Ma_(1.0), u_(1.0), p_(0.0), T_(288.16), rho_(1.225), mu_(1.7894e-5), gamma_(1.4),
              area_(1), length_(1), R_(287.06) {
            Tt_ = T_ * (1 + (gamma_ - 1) / 2.0 * sqr(Ma_));
        }

        referenceValues(const referenceValues &) = delete;
        referenceValues &operator=(const referenceValues &) = delete;

        /**
         * \brief Construct from controller.
         */
        referenceValues(const controller &cont);

        /**
         * \brief Destructor.
         */
        inline ~referenceValues() noexcept;

        // Member functions

        void reset(const controller &cont);

        /**
         * \brief Return reference Mach number.
         */
        hur_nodiscard inline real Ma() const noexcept;

        /**
         * \brief Return reference velocity [m/s].
         */
        hur_nodiscard inline real vMag() const noexcept;

        /**
         * \brief Return reference pressure [Pascal].
         */
        hur_nodiscard inline real p() const noexcept;

        /**
         * \brief Return reference temperature [K].
         */
        hur_nodiscard inline real T() const noexcept;

        /**
         * \brief Return reference stagnation temperature [K].
         */
        hur_nodiscard inline real Tt() const noexcept;

        /**
         * \brief Return reference density [kg/m^3].
         */
        hur_nodiscard inline real rho() const noexcept;

        /**\brief Return reference viscosity [kg/(m s)]*/
        hur_nodiscard inline real mu() const noexcept;

        /**
         * \brief Return ratio of specific heats.
         */
        hur_nodiscard inline real gama() const noexcept;

        /** \brief The reference area, which is used to compute the force and the moment coefficients. [m^2] */
        hur_nodiscard inline real area() const noexcept;

        /** \brief The reference length, which is used to compute the moment coefficient. [m] */
        hur_nodiscard inline real length() const noexcept;
    };

} // namespace OpenHurricane

#include "referenceValues.inl"