/*!
 * \file randoms.hpp
 * \brief Header of randoms.
 *       The subroutines and functions are in the <i>randoms.cpp</i> file.
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
#include "realArray.hpp"
#include <chrono>
#include <random>

namespace OpenHurricane {
    class randoms {
    private:
        std::normal_distribution<real> N_;

    public:
        /**
         * \brief Disallow null constructor.
         */
        randoms() = delete;

        /**
         * \brief Construct from mean and standard deviation.
         */
        randoms(const real mean, const real stddev);

        /**
         * \brief Destructor.
         */
        inline ~randoms() noexcept {}

        /**
         * \brief Operator to generate an array of normally distributed ramdon datas.
         * \param[in] n - The size of array
         * \param[in] seed - The seed for random engine
         */
        hur_nodiscard realArray operator()(
            const integer n,
            const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count());
    };

} // namespace OpenHurricane