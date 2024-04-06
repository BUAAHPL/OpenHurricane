/*!
 * \file sinteger.hpp
 * \brief Header of short int
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
#include <climits>
#include <cstdint>
#include <cstdlib>

#include "HurMPIBase.hpp"
#include "feature.hpp"
#include "string.hpp"
#include "setClassName.hpp"
namespace OpenHurricane {

    using sinteger = short;
    using usinteger = unsigned short;

    /**\brief Template specialization for feature<short_t>.*/
    template <> class feature<sinteger> {
        sinteger p_;

    public:
        using elementType = sinteger;

        /**\brief MPI datatype.*/
        static constexpr HurMPIBase::Datatype MPIType = MPI_SHORT;

        static constexpr int rank = 0;
        static constexpr int nElements_ = 1;

        /*!\brief Construct from primitive.*/
        constexpr inline explicit feature(const elementType &p) : p_(p) {}

        /*!\brief Construct from IStringStream.*/
        feature(IStringStream &is);

        /*!\brief Access to the short value.*/
        constexpr inline operator elementType() const noexcept { return p_; }
        constexpr inline operator elementType &() noexcept { return p_; }
    };

    hur_nodiscard inline sinteger mag(const sinteger l) {
        return (sinteger)std::abs(l);
    }

    template <> hur_nodiscard inline std::string toString(const sinteger &is) {
        return std::to_string(is);
    }

    static const sinteger shortMax = SHRT_MAX;
    static const sinteger shortMin = SHRT_MIN;

    /**\brief Template specialization for feature<short_t>.*/
    template <> class feature<usinteger> {
        usinteger p_;

    public:
        using elementType = usinteger;

        /**\brief MPI datatype.*/
        static constexpr HurMPIBase::Datatype MPIType = MPI_UNSIGNED_SHORT;

        static constexpr int rank = 0;
        static constexpr int nElements_ = 1;

        /*!\brief Construct from primitive.*/
        constexpr inline explicit feature(const usinteger &p) : p_(p) {}

        /*!\brief Construct from IStringStream.*/
        feature(IStringStream &is);

        /*!\brief Access to the short value.*/
        constexpr inline operator usinteger() const noexcept{ return p_; }
        constexpr inline operator usinteger &() noexcept { return p_; }
    };

    hur_nodiscard constexpr inline usinteger mag(const usinteger l) noexcept{
        return l;
    }

    template <> hur_nodiscard inline std::string toString(const usinteger &is) {
        return std::to_string(is);
    }

    static const usinteger ushortMax = USHRT_MAX;

} // namespace OpenHurricane