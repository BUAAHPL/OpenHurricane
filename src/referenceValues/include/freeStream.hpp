/*!
 * \file basicFreeStreams.hpp
 * \brief Header of free streams
 *       The subroutines and functions are in the <i>basicFreeStreams.cpp</i> file.
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

#include "controller.hpp"
#include "dataStructure.hpp"
#include "realArray.hpp"

namespace OpenHurricane {
    class flowModel;

    class freeStream {
    private:
        const flowModel &flow_;

        /** \brief Mach number of free stream. */
        real Ma_;

        /** \brief Velocity of free stream. */
        vector v_;

        /** \brief Magnitude of free stream velocity. */
        real vMag_;

        /** \brief Static pressure of free stream. */
        real p_;

        /** \brief Static temperature of free stream. */
        real T_;

        /** \brief Density of free stream. */
        real rho_;

        real mu_;
        real gamma_;

        realArray yi_;

        void autoSetFreeStream(const controller &bbCont);

        void setFreeStream(const controller &bcCont);

        void getSpeciesMassFractions(const controller &bcCont);

    public:
        /**
         * \brief Disallow null constructor.
         */
        freeStream() = delete;

        freeStream(const controller &cont, const flowModel &flow);

        /**
         * \brief Disallow copy constructor.
         */
        freeStream(const freeStream &) = delete;
        freeStream &operator=(const freeStream &fs) = delete;

        /**
         * \brief Destructor.
         */
        virtual ~freeStream() noexcept {}

        /** \brief Mach number of free stream. */
        hur_nodiscard inline real Ma() const noexcept;

        /** \brief Velocity of free stream. */
        hur_nodiscard inline const vector &v() const noexcept;

        /** \brief Magnitude of free stream velocity. */
        hur_nodiscard inline real vMag() const noexcept;

        /** \brief Static pressure of free stream. */
        hur_nodiscard inline real p() const noexcept;

        /** \brief Static temperature of free stream. */
        hur_nodiscard inline real T() const noexcept;

        /** \brief Density of free stream. */
        hur_nodiscard inline real rho() const noexcept;

        /**\brief Return molecular viscosity*/
        hur_nodiscard inline real mu() const noexcept;

        hur_nodiscard inline real gama() const noexcept;

        hur_nodiscard inline const realArray &yi() const noexcept;
    };

} // namespace OpenHurricane

#include "freeStream.inl"