/*!
 * \file doubleReal.hpp
 * \brief Header of double real
 *       The subroutines and functions are in the <i>doubleReal.cpp</i> file.
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
#include "feature.hpp"
#include "integer.hpp"
#include "string.hpp"
#include <iomanip>
#include <math.h>
#include <stdlib.h>
namespace OpenHurricane {
    using doubleReal = double;

    static constexpr doubleReal doubleRealLarge = 1.0e+15;
    static constexpr doubleReal doubleRealVeryLarge = 1.0e+300;
    static constexpr doubleReal doubleRealRootVeryLarge = 1.0e+150;
    static constexpr doubleReal doubleRealTiny = 1.0e-15;
    static constexpr doubleReal doubleRealRootTiny = 3.0e-8;
    static constexpr doubleReal doubleRealVeryTiny = 1.0e-300;
    static constexpr doubleReal doubleRealRootVeryTiny = 1.0e-150;

    // Read whole of buf as a real. Return true if succesful.
    inline bool readDoubleReal(const char *buf, doubleReal &s) {
        char *endPtr;
        s = strtod(buf, &endPtr);

        return (*endPtr == '\0');
    }

    hur_nodiscard inline doubleReal mag(const doubleReal s) {
        return std::fabs(s);
    }

    /**\brief template specialisation for feature<doubleReal>.*/
    template <> class feature<doubleReal> {
        doubleReal p_;

    public:
        using elementType = doubleReal;
        using integerType = integer;

        /**\brief MPI datatype.*/
        static constexpr HurMPIBase::Datatype MPIType = MPI_DOUBLE;

        /**\brief Tecplot data format.*/
        static constexpr int dataFormat = 2;

        /**\brief Rank of doubleReal is 0.*/
        static constexpr int rank = 0;

        /**\brief Number of components in doubleReal is 1.*/
        static constexpr int nElements_ = 1;

        static constexpr std::streamsize precision = 17;

        /**\brief The Precision used in cgns.*/
        static constexpr int cgnsPrecision = 64;

        static hur_nodiscard constexpr inline bool isFloat() noexcept;
        static hur_nodiscard constexpr inline bool isDouble() noexcept;

        // Constructors

        /*!\brief Construct from primitive.*/
        constexpr inline explicit feature(const doubleReal &p) : p_(p) {}

        /*!\brief Construct from IStringStream.*/
        feature(IStringStream &is);

        // Member Functions

        /*!\brief Access to the doubleReal value.*/
        hur_nodiscard constexpr inline operator doubleReal() const noexcept { return p_; }

        /*!\brief Access to the doubleReal value.*/
        hur_nodiscard constexpr inline operator doubleReal &() noexcept { return p_; }
    };

    hur_nodiscard inline doubleReal fabs(const doubleReal s) {
        return std::fabs(s);
    }

    hur_nodiscard inline doubleReal exp(const doubleReal s) {
        return std::exp(s);
    }

    hur_nodiscard inline doubleReal sin(const doubleReal s) {
        return std::sin(s);
    }
    hur_nodiscard inline doubleReal sinh(const doubleReal s) {
        return std::sinh(s);
    }

    hur_nodiscard inline doubleReal cos(const doubleReal s) {
        return std::cos(s);
    }

    hur_nodiscard inline doubleReal cosh(const doubleReal s) {
        return std::cosh(s);
    }

    hur_nodiscard inline doubleReal tan(const doubleReal s) {
        return std::tan(s);
    }

    hur_nodiscard inline doubleReal tanh(const doubleReal s) {
        return std::tanh(s);
    }

    hur_nodiscard inline doubleReal asin(const doubleReal s) {
        return std::asin(s);
    }

    hur_nodiscard inline doubleReal asinh(const doubleReal s) {
        return std::asinh(s);
    }

    hur_nodiscard inline doubleReal acos(const doubleReal s) {
        return std::acos(s);
    }

    hur_nodiscard inline doubleReal acosh(const doubleReal s) {
        return std::acosh(s);
    }

    hur_nodiscard inline doubleReal atan(const doubleReal s) {
        return std::atan(s);
    }

    hur_nodiscard inline doubleReal atan2(const doubleReal y, const doubleReal x) {
        return std::atan2(y, x);
    }

    hur_nodiscard inline doubleReal atanh(const doubleReal s) {
        return std::atanh(s);
    }

    hur_nodiscard inline doubleReal log(const doubleReal s) {
        return std::log(s);
    }

    hur_nodiscard inline doubleReal log10(const doubleReal s) {
        return std::log10(s);
    }

    hur_nodiscard inline doubleReal sqrt(const doubleReal s) {
        return std::sqrt(s);
    }

    hur_nodiscard inline doubleReal hypot(const doubleReal x, const doubleReal y) {
        return std::hypot(x, y);
    }

    hur_nodiscard inline doubleReal ceil(const doubleReal s) {
        return std::ceil(s);
    }

    hur_nodiscard inline doubleReal floor(const doubleReal s) {
        return std::floor(s);
    }

    hur_nodiscard inline doubleReal fmod(const doubleReal x, const doubleReal y) {
        return std::fmod(x, y);
    }

    inline constexpr doubleReal &setComponent(doubleReal &s, const int) noexcept {
        return s;
    }

    inline constexpr doubleReal component(const doubleReal s, const int) noexcept {
        return s;
    }

    hur_nodiscard constexpr inline doubleReal
    min(const doubleReal &_Left, const doubleReal &_Right) noexcept(noexcept(_Right < _Left)) {
        return std::min(_Left, _Right);
        //return (l < k) ? l : k;
    }

    hur_nodiscard constexpr inline doubleReal
    max(const doubleReal &_Left, const doubleReal &_Right) noexcept(noexcept(_Left < _Right)) {
        return std::max(_Left, _Right);
    }

    // Return 1 if s is positive but not 0
    hur_nodiscard constexpr inline doubleReal notZeroD(const doubleReal s) noexcept {
        return (s != 0) ? s : doubleRealVeryTiny;
    }

    // Return 1 if s is positive or 0 otherwise -1
    hur_nodiscard constexpr inline doubleReal sign(const doubleReal s) noexcept {
        return (s >= 0) ? doubleReal(1) : doubleReal(-1);
    }

    hur_nodiscard inline doubleReal limit(const doubleReal s1, const doubleReal s2) {
        return (mag(s1) < mag(s2)) ? s1 : doubleReal(0.0);
    }

    hur_nodiscard constexpr inline doubleReal magSqr(const doubleReal s) noexcept {
        return s * s;
    }

    // Square
    hur_nodiscard constexpr inline doubleReal sqr(const doubleReal s) noexcept {
        return s * s;
    }

    hur_nodiscard constexpr inline doubleReal pow3(const doubleReal s) noexcept {
        return s * sqr(s);
    }

    hur_nodiscard constexpr inline doubleReal pow4(const doubleReal s) noexcept {
        return sqr(sqr(s));
    }

    hur_nodiscard constexpr inline doubleReal pow5(const doubleReal s) noexcept {
        return s * pow4(s);
    }

    hur_nodiscard constexpr inline doubleReal pow6(const doubleReal s) noexcept {
        return pow3(sqr(s));
    }

    hur_nodiscard inline doubleReal pow025(const doubleReal s) noexcept {
        return sqrt(sqrt(s));
    }

    hur_nodiscard constexpr inline doubleReal inv(const doubleReal s) noexcept {
        return doubleReal(1.0) / s;
    }

    hur_nodiscard constexpr inline doubleReal dot(const doubleReal s1,
                                                  const doubleReal s2) noexcept {
        return s1 * s2;
    }

    hur_nodiscard constexpr inline doubleReal componentAdd(const doubleReal s1,
                                                           const doubleReal s2) noexcept {
        return s1 + s2;
    }

    hur_nodiscard constexpr inline doubleReal componentSubtract(const doubleReal s1,
                                                                const doubleReal s2) noexcept {
        return s1 - s2;
    }

    hur_nodiscard constexpr inline doubleReal componentMultiply(const doubleReal s1,
                                                                const doubleReal s2) noexcept {
        return s1 * s2;
    }

    hur_nodiscard inline doubleReal componentPow(const doubleReal s1, const doubleReal s2) {
        return std::pow(s1, s2);
    }

    hur_nodiscard constexpr inline doubleReal componentDivide(const doubleReal s1,
                                                              const doubleReal s2) noexcept {
        return s1 / s2;
    }

    hur_nodiscard inline constexpr const doubleReal &componentMax(const doubleReal &s) noexcept {
        return s;
    }

    hur_nodiscard constexpr inline const doubleReal componentMax(const doubleReal &s1,
                                                                 const doubleReal &s2) noexcept {
        return max(s1, s2);
    }

    hur_nodiscard inline constexpr const doubleReal &componentMin(const doubleReal &s) noexcept {
        return s;
    }

    hur_nodiscard constexpr inline const doubleReal componentMin(const doubleReal &s1,
                                                                 const doubleReal &s2) noexcept {
        return min(s1, s2);
    }

    hur_nodiscard inline constexpr const doubleReal &componentAv(const doubleReal &s) noexcept {
        return s;
    }

    hur_nodiscard constexpr inline doubleReal componentSqr(const doubleReal &s) noexcept {
        return sqr(s);
    }

    hur_nodiscard inline doubleReal componentMag(const doubleReal &s) {
        return mag(s);
    }

    hur_nodiscard constexpr inline doubleReal componentSign(const doubleReal &s) {
        return sign(s);
    }

    hur_nodiscard inline doubleReal sqrtSumSqr(const doubleReal a, const doubleReal b) {
        doubleReal maga = mag(a);
        doubleReal magb = mag(b);

        if (maga > magb) {
            return maga * sqrt(doubleReal(1.0) + sqr(magb / maga));
        } else {
            return magb < doubleRealVeryTiny ? doubleReal(0.0)
                                             : magb * sqrt(doubleReal(1.0) + sqr(maga / magb));
        }
    }

    // Stabilisation around zero for division
    hur_nodiscard constexpr inline doubleReal stabilise(const doubleReal s,
                                                        const doubleReal small) noexcept {
        if (s >= 0) {
            return s + small;
        } else {
            return s - small;
        }
    }

    hur_forceinline hur_nodiscard doubleReal pow(const doubleReal x, const doubleReal y) {
        return std::pow(x, y);
    }

    hur_forceinline hur_nodiscard doubleReal pow(const doubleReal x, const integer y) {
        return std::pow(x, y);
    }

    template <> hur_nodiscard inline std::string toString(const doubleReal &is) {
        std::stringstream ss;
        ss << std::setprecision(15) << is;
        return ss.str();
    }

    hur_nodiscard constexpr inline bool feature<doubleReal>::isFloat() noexcept {
        return false;
    }

    hur_nodiscard constexpr inline bool feature<doubleReal>::isDouble() noexcept {
        return true;
    }

} // namespace OpenHurricane