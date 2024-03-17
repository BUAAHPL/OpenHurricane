/*!
 * \file clock.hpp
 * \brief Headers of the clock.
 *        The subroutines and functions are in the <i>clock.cpp</i> file.
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

#include "preset.hpp"
#include <chrono>
#include <ctime>
#include <string>

namespace OpenHurricane {

    /**
     * \brief The class of clock
     */
    class clock {
    private:
        /**\brief Names of the days.*/
        static const char *dayNames[];

    public:
        /**\brief Null constructor which stores the start time.*/
        inline clock() {}

        /**!\brief Destructor.*/
        inline ~clock() noexcept {}

        /**\brief Get the current clock time in seconds.*/
        hur_nodiscard static inline time_t getTime() { return time(reinterpret_cast<time_t *>(0)); }

        /**
         *\brief Return the current wall-clock date/time as a string
         *       format according to ISO-8601 (yyyy-mm-ddThh:mm:ss).
         */
        hur_nodiscard static std::string dateTime();
    };

    /**
     * \brief The class of high resolution clock
     */
    class hrClock {
    private:
        /**\brief Start time.*/
        std::chrono::time_point<std::chrono::high_resolution_clock> startTime_;

        /**\brief Time when clockTimeIncrement() was last called.*/
        mutable std::chrono::time_point<std::chrono::high_resolution_clock> lastTime_;

        /**\brief Latest time from either elapsedClockTime() or clockTimeIncrement().*/
        mutable std::chrono::time_point<std::chrono::high_resolution_clock> newTime_;

    public:
        /**\brief Null constructor which stores the start time.*/
        inline hrClock()
            : startTime_(getClockTime()), lastTime_(startTime_), newTime_(startTime_) {}

        /**!\brief Destructor.*/
        inline ~hrClock() noexcept {}

        // Member Functions

        /**\brief Get the current hrClock time in seconds.*/
        hur_nodiscard static inline auto getClockTime()
            -> std::chrono::time_point<std::chrono::high_resolution_clock> {
            return std::chrono::high_resolution_clock::now();
        }

        /**\brief Returns wall-clock time from clock instantiation in seconds.*/
        hur_nodiscard double elapsedClockTime() const;

        /**\brief Returns wall-clock time from last call of clockTimeIncrement() in seconds.*/
        hur_nodiscard double clockTimeIncrement() const;
    };
} // namespace OpenHurricane