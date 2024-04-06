/*!
 * \file expertSystemAdaptCFL.hpp
 * \brief Header of expert-ystem adapt CFL number.
 *       The subroutines and functions are in the <i>expertSystemAdaptCFL.cpp</i> file.
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
#include "CFL.hpp"
#include "cellArrays.hpp"
#include "iteration.hpp"

namespace OpenHurricane {
    /*!\brief The linear progression of CFL number.*/
    class expertSystemAdaptCFL : public CFL {
    private:
        mutable real CFLFactor0_;
        real CFLFactor0Max_;

        real breakdownFactor_;
        real divergenceFactor_;
        real residualJumpThreshold_;
        real updateFactor_;

        integer interval_;
        mutable integer countUnchange_;

        mutable integer countResUp_;

        mutable realArray m2_;
        mutable real lastMaxRes_;
        mutable real lastAveRes_;

        mutable real lastEMaxRes_;
        mutable real lastEAveRes_;

        mutable bool isBreakdown_;

    public:
        declareClassNames;

        /*!\brief Disallow null constructor.*/
        expertSystemAdaptCFL() = delete;

        /*!\brief Construct from controller.*/
        expertSystemAdaptCFL(const iteration &iter, const runtimeMesh &mesh,
                             const controller &cont);

        /*!\brief Disallow copy constructor.*/
        expertSystemAdaptCFL(const expertSystemAdaptCFL &cfl) = delete;

        /*!\brief Destructor.*/
        virtual ~expertSystemAdaptCFL() noexcept {}

        /*!\brief Return the CFL number for current step.
         * \param[in] cstep - The current step.
         * \return The CFL number for current step.
         * \retval A real value.
         */
        virtual real getCFL() const;

        virtual void setRelativeCorrection(const real relCorre) const;
        virtual void setBreakdowns() const;

    private:
        /**
         * \brief If the calculation is broken down, then reduce the CFL number.
         */
        void breakdowns() const;

        /**
         * \brief If the calculation is divergence, then reduce the CFL number.
         */
        void divergence() const;

        /**
         * \brief If the convergence speed is slow, then increase the CFL number.
         */
        void slowConvergence() const;

        /*!\brief Print CFL number */
        void printCFL(const real cfl) const;

        /**
         * \brief Should activate adapting CFL based on the continuous convergence redisuals steps.
         */
        bool activateAdaptCFL(const integer icount) const;
    };
} // namespace OpenHurricane
