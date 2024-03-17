/*!
 * \file bools.hpp
 * \brief Header of bool
 *       The subroutines and functions are in the <i>bools.cpp</i> file.
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
#include "string.hpp"
#include "setClassName.hpp"
namespace OpenHurricane {

    /**\brief Template specialization for feature<bool>.*/
    template <> class feature<bool> {
    private:
        bool p_;

    public:
        using elementType = bool;
        static constexpr int rank = 0;
        static constexpr int nElements_ = 1;
#ifdef MPI_CXX_BOOL
        // MPI datatype
        static constexpr HurMPIBase::Datatype MPIType = MPI_CXX_BOOL;
        //static const HurMPIBase::Datatype MPIType = MPI_CXX_BOOL;
#else
        // MPI datatype
        static constexpr HurMPIBase::Datatype MPIType = MPI_C_BOOL;
        //static const HurMPIBase::Datatype MPIType = MPI_C_BOOL;
#endif // MPI_CXX_BOOL

        /*!\brief Construct from primitive.*/
        constexpr inline explicit feature(const bool &p) : p_(p) {}

        /*!\brief Construct from IStringStream.*/
        feature(IStringStream &is);

        /*!\brief Access to the bool value.*/
        constexpr inline operator bool() const noexcept { return p_; }
        constexpr inline operator bool &() noexcept { return p_; }
    };

    template <> hur_nodiscard inline std::string toString(const bool &is) {
        std::string s = "False";
        if (is) {
            s = "True";
        }
        return s;
    }

} // namespace OpenHurricane
