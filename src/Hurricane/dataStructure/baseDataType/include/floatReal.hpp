/*!
 * \file floatReal.hpp
 * \brief Header of float real
 *       The subroutines and functions are in the <i>floatReal.cpp</i> file.
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
#include "feature.hpp"
#include "integer.hpp"
#include "string.hpp"
#include <cmath>
#include <stdlib.h>

namespace OpenHurricane {
    using floatReal = float;

    static constexpr floatReal floatRealLarge = 1.0e+6f;
    static constexpr floatReal floatRealVeryLarge = 1.0e+37f;
    static constexpr floatReal floatRealRootVeryLarge = 1.0e+18f;
    static constexpr floatReal floatRealTiny = 1.0e-6f;
    static constexpr floatReal floatRealRootTiny = 3.0e-3f;
    static constexpr floatReal floatRealVeryTiny = 1.0e-37f;
    static constexpr floatReal floatRealRootVeryTiny = 1.0e-18f;

    // Read whole of buf as a real. Return true if succesful.
    inline bool readFloatReal(const char *buf, floatReal &s) {
        char *endPtr;
        s = strtof(buf, &endPtr);

        return (*endPtr == '\0');
    }

    /**\brief template specialisation for feature<floatReal>.*/
    template <> class feature<floatReal> {
        floatReal p_;

    public:
        using elementType = floatReal;

        using integerType = integer;

        /**\brief MPI datatype.*/
        static constexpr HurMPIBase::Datatype MPIType = MPI_FLOAT;

        /**\brief Tecplot data format.*/
        static constexpr int dataFormat = 1;

        /**\brief Rank of floatReal is 0.*/
        static constexpr int rank = 0;

        /**\brief Number of components in floatReal is 1.*/
        static constexpr int nElements_ = 1;

        static constexpr std::streamsize precision = 6;

        /**\brief The Precision used in cgns.*/
        static constexpr int cgnsPrecision = 32;

        static hur_nodiscard constexpr inline bool isFloat() noexcept;
        static hur_nodiscard constexpr inline bool isDouble() noexcept;

        // Constructors

        /*!\brief Construct from primitive.*/
        constexpr inline explicit feature(const floatReal &p) : p_(p) {}

        /*!\brief Construct from IStringStream.*/
        feature(IStringStream &is);

        // Member Functions

        /*!\brief Access to the floatReal value.*/
        hur_nodiscard constexpr inline operator floatReal() const noexcept { return p_; }

        /*!\brief Access to the floatReal value.*/
        hur_nodiscard constexpr inline operator floatReal &() noexcept { return p_; }
    };

    hur_nodiscard inline floatReal fabs(const floatReal s) {
        return std::fabs(s);
    }

    hur_nodiscard inline floatReal exp(const floatReal s) {
        return std::exp(s);
    }

    hur_nodiscard inline floatReal sin(const floatReal s) {
        return std::sin(s);
    }
    hur_nodiscard inline floatReal sinh(const floatReal s) {
        return std::sinh(s);
    }

    hur_nodiscard inline floatReal cos(const floatReal s) {
        return std::cos(s);
    }

    hur_nodiscard inline floatReal cosh(const floatReal s) {
        return std::cosh(s);
    }

    hur_nodiscard inline floatReal tan(const floatReal s) {
        return std::tan(s);
    }

    hur_nodiscard inline floatReal tanh(const floatReal s) {
        return std::tanh(s);
    }

    hur_nodiscard inline floatReal asin(const floatReal s) {
        return std::asin(s);
    }

    hur_nodiscard inline floatReal asinh(const floatReal s) {
        return std::asinh(s);
    }

    hur_nodiscard inline floatReal acos(const floatReal s) {
        return std::acos(s);
    }

    hur_nodiscard inline floatReal acosh(const floatReal s) {
        return std::acosh(s);
    }

    hur_nodiscard inline floatReal atan(const floatReal s) {
        return std::atan(s);
    }

    hur_nodiscard inline floatReal atan2(const floatReal y, const floatReal x) {
        return std::atan2(y, x);
    }

    hur_nodiscard inline floatReal atanh(const floatReal s) {
        return std::atanh(s);
    }

    hur_nodiscard inline floatReal log(const floatReal s) {
        return std::log(s);
    }

    hur_nodiscard inline floatReal log10(const floatReal s) {
        return std::log10(s);
    }

    hur_nodiscard inline floatReal sqrt(const floatReal s) {
        return std::sqrt(s);
    }

    hur_nodiscard inline floatReal hypot(const floatReal x, const floatReal y) {
        return std::hypot(x, y);
    }

    hur_nodiscard inline floatReal ceil(const floatReal s) {
        return std::ceil(s);
    }

    hur_nodiscard inline floatReal floor(const floatReal s) {
        return std::floor(s);
    }

    hur_nodiscard inline floatReal mag(const floatReal s) {
        return fabsf(s);
    }

    hur_nodiscard inline floatReal fmod(const floatReal x, const floatReal y) {
        return std::fmod(x, y);
    }

    inline constexpr floatReal &setComponent(floatReal &s, const int) noexcept {
        return s;
    }

    inline constexpr floatReal component(const floatReal s, const int) noexcept {
        return s;
    }

    hur_nodiscard constexpr inline floatReal
    min(const floatReal &_Left, const floatReal &_Right) noexcept(noexcept(_Right < _Left)) {
        return std::min(_Left, _Right);
    }

    hur_nodiscard constexpr inline floatReal
    max(const floatReal &_Left, const floatReal &_Right) noexcept(noexcept(_Left < _Right)) {
        return std::max(_Left, _Right);
    }

    // Return 1 if s is positive but not 0
    hur_nodiscard constexpr inline floatReal notZeroD(const floatReal s) noexcept {
        return (s != 0) ? s : floatRealVeryTiny;
    }

    // Return 1 if s is positive or 0 otherwise -1
    hur_nodiscard constexpr inline floatReal sign(const floatReal s) noexcept {
        return (s >= 0) ? floatReal(1) : floatReal(-1);
    }

    hur_nodiscard inline floatReal limit(const floatReal s1, const floatReal s2) {
        return (mag(s1) < mag(s2)) ? s1 : floatReal(0.0);
    }

    hur_nodiscard constexpr inline floatReal magSqr(const floatReal s) noexcept {
        return s * s;
    }

    // Square
    hur_nodiscard constexpr inline floatReal sqr(const floatReal s) noexcept {
        return s * s;
    }

    hur_nodiscard constexpr inline floatReal pow3(const floatReal s) noexcept {
        return s * sqr(s);
    }

    hur_nodiscard constexpr inline floatReal pow4(const floatReal s) noexcept {
        return sqr(sqr(s));
    }

    hur_nodiscard constexpr inline floatReal pow5(const floatReal s) noexcept {
        return s * pow4(s);
    }

    hur_nodiscard constexpr inline floatReal pow6(const floatReal s) noexcept {
        return pow3(sqr(s));
    }

    hur_nodiscard inline floatReal pow025(const floatReal s) noexcept {
        return sqrt(sqrt(s));
    }

    hur_nodiscard inline floatReal inv(const floatReal s) noexcept {
        return floatReal(1.0) / s;
    }

    hur_nodiscard constexpr inline floatReal dot(const floatReal s1, const floatReal s2) noexcept {
        return s1 * s2;
    }

    hur_nodiscard constexpr inline floatReal componentAdd(const floatReal s1,
                                                          const floatReal s2) noexcept {
        return s1 + s2;
    }

    hur_nodiscard constexpr inline floatReal componentSubtract(const floatReal s1,
                                                               const floatReal s2) noexcept {
        return s1 - s2;
    }

    hur_nodiscard constexpr inline floatReal componentMultiply(const floatReal s1,
                                                               const floatReal s2) noexcept {
        return s1 * s2;
    }

    hur_nodiscard inline floatReal componentPow(const floatReal s1, const floatReal s2) {
        return std::pow(s1, s2);
    }

    hur_nodiscard constexpr inline floatReal componentDivide(const floatReal s1,
                                                             const floatReal s2) noexcept {
        return s1 / s2;
    }

    hur_nodiscard constexpr inline const floatReal &componentMax(const floatReal &s) noexcept {
        return s;
    }

    hur_nodiscard constexpr inline const floatReal componentMax(const floatReal &s1,
                                                                const floatReal &s2) {
        return max(s1, s2);
    }

    hur_nodiscard constexpr inline const floatReal &componentMin(const floatReal &s) noexcept {
        return s;
    }

    hur_nodiscard constexpr inline const floatReal componentMin(const floatReal &s1,
                                                                const floatReal &s2) {
        return min(s1, s2);
    }

    hur_nodiscard inline constexpr const floatReal &componentAv(const floatReal &s) noexcept {
        return s;
    }

    hur_nodiscard constexpr inline floatReal componentSqr(const floatReal &s) noexcept {
        return sqr(s);
    }

    hur_nodiscard inline floatReal componentMag(const floatReal &s) {
        return mag(s);
    }

    hur_nodiscard constexpr inline floatReal componentSign(const floatReal &s) {
        return sign(s);
    }

    hur_nodiscard inline floatReal sqrtSumSqr(const floatReal a, const floatReal b) {
        floatReal maga = mag(a);
        floatReal magb = mag(b);

        if (maga > magb) {
            return maga * sqrt(floatReal(1.0) + sqr(magb / maga));
        } else {
            return magb < floatRealVeryTiny ? floatReal(0.0)
                                            : magb * sqrt(floatReal(1.0) + sqr(maga / magb));
        }
    }

    // Stabilisation around zero for division
    hur_nodiscard constexpr inline floatReal stabilise(const floatReal s, const floatReal small) {
        if (s >= 0) {
            return s + small;
        } else {
            return s - small;
        }
    }
} // namespace OpenHurricane

namespace OpenHurricane {

    hur_nodiscard inline floatReal pow(const floatReal x, const floatReal y) {
        //return std::powf(x, y);
        return std::pow(x, y);
    }

    hur_nodiscard inline floatReal pow(const floatReal x, const integer y) {
        return (floatReal)std::pow(x, y);
    }

    template <> hur_nodiscard inline std::string toString(const floatReal &is) {
        return std::to_string(is);
    }

    hur_nodiscard constexpr inline bool feature<floatReal>::isFloat() noexcept {
        return true;
    }

    hur_nodiscard constexpr inline bool feature<floatReal>::isDouble() noexcept {
        return false;
    }

} // namespace OpenHurricane
