/*!
 * \file int64.hpp
 * \brief Header of int64
 *       The subroutines and functions are in the <i>int64.cpp</i> file.
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
#include <climits>
#include <cstdint>
#include <cstdlib>

#include "feature.hpp"
#include "setClassName.hpp"
#include "string.hpp"
namespace OpenHurricane {
    /**
     * \brief Template specialization for feature<int64_t>.
     */
    template <> class feature<int64_t> {
        int64_t p_;

    public:
        using elementType = int64_t;

        /**
         * \brief MPI datatype.
         */
        static constexpr HurMPIBase::Datatype MPIType = MPI_INT64_T;

        static constexpr int rank = 0;
        static constexpr int nElements_ = 1;

        /**\brief data format.*/
        static constexpr int dataFormat = 3;

        static constexpr int64_t min = INT64_MIN;
        static constexpr int64_t max = INT64_MAX;
        static constexpr int64_t rootMax = INT64_MIN;
        static constexpr int64_t rootMin = INT64_MAX;

        static constexpr std::streamsize precision = 6;

        /*!\brief Construct from primitive.*/
        constexpr inline explicit feature(const int64_t &p) : p_(p) {}

        /*!\brief Construct from IStringStream.*/
        feature(IStringStream &is);

        /*!\brief Access to the int64_t value.*/
        constexpr inline operator int64_t() const noexcept { return p_; }
        constexpr inline operator int64_t &() noexcept { return p_; }
    };

    hur_nodiscard inline int64_t mag(const int64_t l) {
        return std::abs(l);
    }

    template <> hur_nodiscard inline std::string toString(const int64_t &is) {
        return std::to_string(is);
    }

    /**
     * \brief Template specialization for feature<uint64_t>.
     */
    template <> class feature<uint64_t> {
        uint64_t p_;

    public:
        using elementType = uint64_t;

        /**
         * \brief MPI datatype..
         */
        static constexpr HurMPIBase::Datatype MPIType = MPI_UINT64_T;

        static constexpr int rank = 0;
        static constexpr int nElements_ = 1;

        /**\brief data format.*/
        static constexpr int dataFormat = 3;

        static constexpr uint64_t min = 0;
        static constexpr uint64_t max = UINT64_MAX;
        static constexpr uint64_t rootMax = 0;
        static constexpr uint64_t rootMin = UINT64_MAX;

        static constexpr std::streamsize precision = 6;

        /*!\brief Construct from primitive.*/
        constexpr inline explicit feature(const uint64_t &p) : p_(p) {}

        /*!\brief Construct from IStringStream.*/
        feature(IStringStream &is);

        /*!\brief Access to the int64_t value.*/
        constexpr inline operator uint64_t() const noexcept { return p_; }
        constexpr inline operator uint64_t &() noexcept { return p_; }
    };

    hur_nodiscard constexpr inline uint64_t mag(const uint64_t l) noexcept {
        return l;
    }

    template <> hur_nodiscard inline std::string toString(const uint64_t &is) {
        return std::to_string(is);
    }
} // namespace OpenHurricane
