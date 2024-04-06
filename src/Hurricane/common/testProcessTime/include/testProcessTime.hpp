/*!
 * \file testProcessTime.hpp
 * \brief Headers of the testing the process time of OpenHurricane.
 *        The subroutines and functions are in the <i>testProcessTime.cpp</i> file.
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
#include "clock.hpp"
#include "controller.hpp"
#include "fileOsstream.hpp"

namespace OpenHurricane {
    class iteration;

    /**
     * \brief The class of testing the process time of OpenHurricane
     */
    class testProcessTime {
    private:
        /**
         * \brief The file stream for output time info.
         */
        mutable fileOsstream fosTestPT_;

        /**
         * \brief The pointer of high resolution clock.
         */
        mutable uniquePtr<hrClock> hClockPtr_;

        mutable std::string titleBuffStr_;
        mutable std::string firstLineBuffStr_;

        mutable integer count_;

    public:
        // Constructors

        /**\brief Disallow null constructor.*/
        testProcessTime() = delete;

        testProcessTime(const iteration &iter, const controller &cont);

        testProcessTime(const iteration &iter, const fileName &fn);

        /**
         * \brief Disallow copy constructor.
         */
        testProcessTime(const testProcessTime &) = delete;

        /**!\brief Destructor.*/
        inline ~testProcessTime() noexcept;

        /**
         * \brief Start timing.
         */
        inline void start() const;

        /**
         * \brief Start timing.
         * \param[in] istep - The step written to file
         */
        void start(const integer istep) const;

        /**
         * \brief Compute wall-clock time from last call of clockTimeIncrement() in seconds.
         */
        hur_nodiscard real clockTimeIncrement() const;

        /**
         * \brief Compute wall-clock time from last call of clockTimeIncrement() in seconds.
         * \param[in] name - The name of this call of clockTimeIncrement (no space)
         */
        void clockTimeIncrement(const char *name) const;

        /**
         * \brief Compute wall-clock time from last call of clockTimeIncrement() in seconds.
         * \param[in] name - The name of this call of clockTimeIncrement (no space)
         */
        void clockTimeIncrement(const char *name, const real time) const;

        /**
         * \brief Stop timing.
         */
        void stop() const;

        /**
         * \brief To clear the pointer of high resolution clock and close the file stream.
         */
        inline void clear() noexcept;

        /**
         * \brief To clear the pointer of high resolution clock.
         */
        inline void resetHZClock() const noexcept;

        /**
         * \brief Disallow bitwise assignment.
         */
        void operator=(const testProcessTime &) = delete;
    };
} // namespace OpenHurricane

#include "testProcessTime.inl"