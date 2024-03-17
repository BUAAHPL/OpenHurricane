/*!
 * \file subIteration.hpp
 * \brief Header of CFD time advance sub-iteration
 *       The subroutines and functions are in the <i>subIteration.cpp</i> file.
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

namespace OpenHurricane {
    /*!\brief The class of sub-iteration for dual-time stepping.*/
    class subIteration {
    private:
        // Private data

        /*!\brief Tha max steps for sub-iteration.*/
        integer maxSubStep_;

        /*!\brief Tha max steps for sub-iteration.*/
        integer minSubStep_;

        /*!\brief Current sub iteration step for dual-time step advance.*/
        integer cSubStep_;

        /*!\brief Tha total steps for entire iteration.*/
        integer totalStep_;

        mutable bool isResConvergence_;

    public:
        subIteration() = default;
        subIteration(const subIteration &) = default;
        subIteration &operator=(const subIteration &) = default;

        /*!\brief Construct from components.*/
        subIteration(const integer maxSubStep, const integer minSubStep = 0);

        /*!\brief Destructor.*/
        inline ~subIteration() noexcept {}

        // Member Functions

        /*!\brief Tha max steps for sub-iteration.*/
        hur_nodiscard inline integer maxSubStep() const noexcept { return maxSubStep_; }

        /*!\brief Tha max steps for sub-iteration.*/
        hur_nodiscard inline integer minSubStep() const noexcept { return minSubStep_; }

        /*!\brief If the sub-iteration is ended.*/
        hur_nodiscard inline bool isLooping() const noexcept { return maxSubStep_ > cSubStep_; }

        /*!\brief If the sub-iteration keeps running.*/
        hur_nodiscard inline bool checkMin() const noexcept { return cSubStep_ <= minSubStep_; }

        /*!\brief Current sub iteration step for dual-time step advance.*/
        hur_nodiscard inline integer cSubStep() const noexcept { return cSubStep_; }

        /*!\brief Total iteration step for dual-time step advance.*/
        hur_nodiscard inline integer totalStep() const noexcept { return totalStep_; }

        /*!\brief Set isResConvergence_ true.*/
        inline void setResConveg() const noexcept { isResConvergence_ = true; }

        /*|!brief Reset the sub-iteration step to zero.*/
        inline void reset() noexcept {
            isResConvergence_ = false;
            cSubStep_ = 0;
        }

        /*|!brief Reset the sub-iteration maxSubstep.*/
        inline void resetMaxStep(const integer mss) noexcept { maxSubStep_ = mss; }

        /*!\brief Sub-iteration.*/
        hur_nodiscard inline bool subIterating() noexcept {
            bool running = (isLooping() && !isResConvergence_) || checkMin();
            if (running) {
                operator++();
            }
            return running;
        }

        /**
         * \brief Is the sub-iteration reaching convergence by checking residuals.
         */
        hur_nodiscard inline bool isResConvergence() const noexcept { return isResConvergence_; }

        inline subIteration &operator++() noexcept {
            totalStep_++;
            cSubStep_++;
            return *this;
        }

        inline subIteration &operator++(int) noexcept { return operator++(); }
    };
} // namespace OpenHurricane
