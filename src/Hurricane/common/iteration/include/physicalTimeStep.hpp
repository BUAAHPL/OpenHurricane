/*!
 * \file physicalTimeStep.hpp
 * \brief Header of physical time step.
 *       The subroutines and functions are in the <i>physicalTimeStep.cpp</i> file.
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
    /*!\brief The class of sub-iteration for dual-time stepping.*/
    class physicalTimeStep {
    private:
        // Private data

        /*!\brief Physical time step [s].*/
        real pTimeStep_;

        /*!\brief Maximum physical time [s].*/
        real maxTime_;

        /*!\brief The total physical time that has been marched [s].*/
        real totalTime_;

        /** \brief If using dynamic time-step. Default is not using.*/
        bool isDynamicTimeStep_;

        bool isDynamicCFLSet_;

        real dyCFL_;

        /**
         * \brief The previous time step [s].
         * lastTimeStep_[0] = dt^n i.e. current time step size.
         * lastTimeStep_[1] = dt^(n-1) i.e. previous time step size.
         * lastTimeStep_[2] = dt^(n-2) i.e.  previous time step size.
         */
        realArray lastTimeStep_;

        void updateLastTimeStep() noexcept;

    public:
        physicalTimeStep() = default;
        physicalTimeStep(const physicalTimeStep &) = default;
        physicalTimeStep &operator=(const physicalTimeStep &) = default;

        /*!\brief Construct from components.*/
        physicalTimeStep(const real timeStep, const real maxTime);

        /*!\brief Construct from components.*/
        physicalTimeStep(const real timeStep, const real maxTime, const real totalTime);

        /*!\brief Destructor.*/
        inline ~physicalTimeStep() noexcept {}

        /*!\brief The current time (has been marching) [s].*/
        hur_nodiscard inline real totalTime() const noexcept { return totalTime_; }

        /*!\brief Physical time step [s].*/
        hur_nodiscard inline real pTimeStep() const noexcept { return pTimeStep_; }

        /*!\brief Maximum physical time [s].*/
        hur_nodiscard inline real maxTime() const noexcept { return maxTime_; }

        /*!\brief If the physical time marches to the max time steps.*/
        hur_nodiscard inline bool end() const noexcept { return totalTime_ >= maxTime_; }

        /**  \brief If using dynamic time-step. Default is not using.*/
        hur_nodiscard inline bool isDynamicTimeStep() const noexcept { return isDynamicTimeStep_; }

        hur_nodiscard inline bool isDynamicCFLSet() const noexcept { return isDynamicCFLSet_; }

        hur_nodiscard inline real dyCFL() const noexcept { return dyCFL_; }

        /** \brief Setting dynamic time-step.*/
        inline bool setDynamicTimeStep() noexcept {
            auto tmp = isDynamicTimeStep_;
            isDynamicTimeStep_ = true;
            return tmp;
        }

        /** \brief Unsetting dynamic time-step.*/
        inline bool unsetDynamicTimeStep() noexcept {
            auto tmp = isDynamicTimeStep_;
            isDynamicTimeStep_ = false;
            isDynamicCFLSet_ = false;
            return tmp;
        }

        inline void setDynamicCFL(const real dycfl) noexcept {
            isDynamicCFLSet_ = true;
            dyCFL_ = dycfl;
        }

        /*\!brief Reset the max time steps [s].*/
        inline void setMaxTime(const real mT) noexcept { maxTime_ = mT; }

        /*!\brief Reset the physical time step [s].*/
        void setTimeStep(const real tS) noexcept;

        /*\!brief Set the current time (has been marching) [s].*/
        inline void setTotalTime(const real tT) noexcept { totalTime_ = tT; }

        hur_nodiscard bool timeMarching() noexcept;

        /**
         * \brief The previous time step [s].
         * lastTimeStep_[0] = dt^n i.e. current time step size.
         * lastTimeStep_[1] = dt^(n-1) i.e. previous time step size.
         * lastTimeStep_[2] = dt^(n-2) i.e.  previous time step size.
         */
        hur_nodiscard inline const realArray &lastTimeStep() const noexcept {
            return lastTimeStep_;
        }

        /**
         * \brief The previous time step [s].
         * lastTimeStep_[0] = dt^n i.e. current time step size.
         * lastTimeStep_[1] = dt^(n-1) i.e. previous time step size.
         * lastTimeStep_[2] = dt^(n-2) i.e.  previous time step size.
         */
        hur_nodiscard inline realArray &lastTimeStep() noexcept { return lastTimeStep_; }

        physicalTimeStep &operator++() noexcept;
        inline physicalTimeStep &operator++(int) noexcept { return operator++(); }
    };
} // namespace OpenHurricane
