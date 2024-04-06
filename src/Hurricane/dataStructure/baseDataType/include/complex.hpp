/*!
 * \file complex.hpp
 * \brief Header of complex
 *       The subroutines and functions are in the <i>complex.cpp</i> file.
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

#include "basicFunctions.hpp"
#include "bools.hpp"
#include "real.hpp"
#include "string.hpp"
#include <complex>

namespace OpenHurricane {

    // Forward declaration of friend functions and operators

    using complex = std::complex<real>;

    hur_nodiscard constexpr inline real magSqr(const complex &c) noexcept {
        return (c.real() * c.real() + c.imag() * c.imag());
    }

    hur_nodiscard inline complex sqr(const complex &c) noexcept {
        return c * c;
    }

    hur_nodiscard inline real mag(const complex &c) {
        return sqrt(magSqr(c));
    }

    hur_nodiscard inline const complex &max(const complex &c1, const complex &c2) {
        if (mag(c1) > mag(c2)) {
            return c1;
        } else {
            return c2;
        }
    }

    hur_nodiscard inline const complex &min(const complex &c1, const complex &c2) {
        if (mag(c1) < mag(c2)) {
            return c1;
        } else {
            return c2;
        }
    }

    hur_nodiscard inline complex limit(const complex &c1, const complex &c2) {
        return complex(limit(c1.real(), c2.real()), limit(c1.imag(), c2.imag()));
    }

    hur_nodiscard constexpr inline complex operator+(const complex &c1,
                                                     const complex &c2) noexcept {
        return complex(c1.real() + c2.real(), c1.imag() + c2.imag());
    }

    hur_nodiscard constexpr inline complex operator-(const complex &c) noexcept {
        return complex(-c.real(), -c.imag());
    }

    hur_nodiscard constexpr inline complex operator-(const complex &c1,
                                                     const complex &c2) noexcept {
        return complex(c1.real() - c2.real(), c1.imag() - c2.imag());
    }

    hur_nodiscard constexpr inline complex operator*(const complex &c1,
                                                     const complex &c2) noexcept {
        return complex(c1.real() * c2.real() - c1.imag() * c2.imag(),
                       c1.imag() * c2.real() + c1.real() * c2.imag());
    }

    hur_nodiscard constexpr inline complex operator/(const complex &c1,
                                                     const complex &c2) noexcept {
        real sqrC2 = magSqr(c2);

        return complex((c1.real() * c2.real() + c1.imag() * c2.imag()) / sqrC2,
                       (c1.imag() * c2.real() - c1.real() * c2.imag()) / sqrC2);
    }

    hur_nodiscard constexpr inline complex operator*(const real s, const complex &c) noexcept {
        return complex(s * c.real(), s * c.imag());
    }

    hur_nodiscard constexpr inline complex operator*(const complex &c, const real s) noexcept {
        return complex(s * c.real(), s * c.imag());
    }

    hur_nodiscard constexpr inline complex operator/(const complex &c, const real s) noexcept {
        return complex(c.real() / s, c.imag() / s);
    }

    hur_nodiscard constexpr inline complex operator/(const real s, const complex &c) noexcept {
        return complex(s / c.real(), s / c.imag());
    }

    template <class Type>
    hur_nodiscard inline std::complex<Type> conjugate(const std::complex<Type> &z) noexcept {
        return std::conj(z);
    }

    /**
     * \brief Template specialization for feature<int64_t>.
     */
    template <> class feature<complex> {
        complex p_;

    public:
        using elementType = real;

        static constexpr HurMPIBase::Datatype MPIType = feature<real>::MPIType;
        static constexpr int nElements_ = 2;
        static constexpr std::streamsize precision = feature<real>::precision;

        /*!\brief Construct from primitive.*/
        constexpr inline explicit feature(const complex &p) : p_(p) {}

        /*!\brief Construct from IStringStream.*/
        feature(IStringStream &is);

        constexpr inline operator const complex &() const noexcept { return p_; }
        constexpr inline operator complex &() noexcept { return p_; }
    };

    template <> hur_nodiscard inline std::string toString(const complex &is) {
        std::string s;
        s = std::to_string(is.real());
        if (is.imag() != real(0)) {
            s += "+" + std::to_string(is.imag()) + "i";
        }
        return s;
    }

} // namespace OpenHurricane
