/*!
 * \file int32.hpp
 * \brief Header of int32
 *       The subroutines and functions are in the <i>int32.cpp</i> file.
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

#include "HurMPIBase.hpp"
#include "feature.hpp"
#include "setClassName.hpp"
#include "string.hpp"
namespace OpenHurricane {
    /**
     * \brief Template specialization for feature<int32_t>.
     */
    template <> class feature<int32_t> {
        int32_t p_;

    public:
        using elementType = int32_t;

        /**
         * \brief MPI datatype.
         */
        static constexpr HurMPIBase::Datatype MPIType = MPI_INT32_T;

        static constexpr int rank = 0;
        static constexpr int nElements_ = 1;

        /**\brief data format.*/
        static constexpr int dataFormat = 4;

        static constexpr int32_t min = INT32_MIN;
        static constexpr int32_t max = INT32_MAX;
        static constexpr int32_t rootMax = INT32_MIN;
        static constexpr int32_t rootMin = INT32_MAX;

        static constexpr std::streamsize precision = 6;

        /*!\brief Construct from primitive.*/
        constexpr inline explicit feature(const int32_t &p) : p_(p) {}

        /*!\brief Construct from IStringStream.*/
        feature(IStringStream &is);

        /*!\brief Access to the int32_t value.*/
        constexpr inline operator int32_t() const noexcept { return p_; }
        constexpr inline operator int32_t &() noexcept { return p_; }
    };

    hur_nodiscard inline int32_t mag(const int32_t l) {
        return std::abs(l);
    }

    template <> hur_nodiscard inline std::string toString(const int32_t &is) {
        return std::to_string(is);
    }

    /**
     * \brief Template specialization for feature<int32_t>.
     */
    template <> class feature<uint32_t> {
        uint32_t p_;

    public:
        using elementType = uint32_t;

        /**
         * \brief MPI datatype.
         */
        static constexpr HurMPIBase::Datatype MPIType = MPI_UINT32_T;

        static constexpr int rank = 0;
        static constexpr int nElements_ = 1;

        /**\brief data format.*/
        static constexpr int dataFormat = 4;

        static constexpr uint32_t min = 0;
        static constexpr uint32_t max = UINT32_MAX;
        static constexpr uint32_t rootMax = 0;
        static constexpr uint32_t rootMin = UINT32_MAX;

        static constexpr std::streamsize precision = 6;

        /*!\brief Construct from primitive.*/
        constexpr inline explicit feature(const uint32_t &p) : p_(p) {}

        /*!\brief Construct from IStringStream.*/
        feature(IStringStream &is);

        /*!\brief Access to the int32_t value.*/
        constexpr inline operator uint32_t() const noexcept { return p_; }
        constexpr inline operator uint32_t &() noexcept { return p_; }
    };

    hur_nodiscard constexpr inline uint32_t mag(const uint32_t l) noexcept {
        return l;
    }

    template <> hur_nodiscard inline std::string toString(const uint32_t &is) {
        return std::to_string(is);
    }
} // namespace OpenHurricane