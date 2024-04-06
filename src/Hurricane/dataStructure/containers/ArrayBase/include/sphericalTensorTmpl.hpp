/*!
 * \file sphericalTensorTmpl.hpp
 * \brief Header of Spherical Tensor
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

#include "vectorTmpl.hpp"

namespace OpenHurricane {
    /**
     * \brief The template class of SphericalTensor.
     */
    template <class Type> class SphericalTensor {
    public:
        using value_type = Type;
        using elementType = value_type;

    private:
        value_type ii_;

    public:
        /**\brief Rank of SphericalTensor is 2*/
        static constexpr int rank = 2;

        static constexpr int nElements_ = 1;
        static const std::streamsize precision;

        static const SphericalTensor I;
        static const SphericalTensor oneThirdI;
        static const SphericalTensor twoThirdsI;

    public:
        constexpr inline SphericalTensor() {}

        constexpr inline SphericalTensor(const OpenHurricane::zero) : ii_(Zero) {}
        constexpr inline SphericalTensor &operator=(const OpenHurricane::zero) {
            ii_ = Zero;
            return *this;
        }

        constexpr inline SphericalTensor(const SphericalTensor &other) : ii_(other.ii_) {}
        constexpr inline SphericalTensor &operator=(const SphericalTensor &other) {
            if (this != std::addressof(other)) {
                ii_ = other.ii_;
            }
            return *this;
        }

        constexpr inline SphericalTensor(const value_type &ii) : ii_(ii) {}

        template <class otherType>
        constexpr inline SphericalTensor(const SphericalTensor<otherType> &other)
            : ii_(static_cast<const value_type &>(other.ii_)) {}

        inline ~SphericalTensor() noexcept {}

        hur_nodiscard constexpr inline const Type &ii() const noexcept { return ii_; }
        hur_nodiscard constexpr inline Type &ii() noexcept { return ii_; }

        /**\brief Transpose*/
        hur_nodiscard constexpr inline const SphericalTensor &transpose() const noexcept {
            return *this;
        }

        hur_nodiscard constexpr inline const Type &operator[](const int i) const noexcept {
            return ii_;
        }
        hur_nodiscard constexpr inline Type &operator[](const int i) noexcept { return ii_; }
    };

    /**\brief Product between a real and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SphericalTensor<Type>
    operator*(const Type s, const SphericalTensor<Type> &st2) {
        return SphericalTensor<Type>(s * st2.ii());
    }

    /**\brief Product between a spherical tensor and a real*/
    template <class Type>
    hur_nodiscard constexpr inline SphericalTensor<Type> operator*(const SphericalTensor<Type> &st2,
                                                                   const Type s) {
        return SphericalTensor<Type>(s * st2.ii());
    }

    /**\brief Inner-product between two spherical tensors*/
    template <class Type>
    hur_nodiscard constexpr inline SphericalTensor<Type>
    operator*(const SphericalTensor<Type> &st1, const SphericalTensor<Type> &st2) {
        return SphericalTensor<Type>(st1.ii() * st2.ii());
    }

    /**\brief Inner-product between two spherical tensors*/
    template <class Type>
    hur_nodiscard constexpr inline SphericalTensor<Type>
    operator&(const SphericalTensor<Type> &st1, const SphericalTensor<Type> &st2) {
        return SphericalTensor<Type>(st1.ii() * st2.ii());
    }

    /**\brief Inner-product between a spherical tensor and a vector*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const SphericalTensor<Type> &st,
                                                          const Vector<Type> &v) {
        return Vector<Type>(st.ii() * v.x(), st.ii() * v.y(), st.ii() * v.z());
    }

    /**\brief Inner-product between a spherical tensor and a vector*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator&(const SphericalTensor<Type> &st,
                                                          const Vector<Type> &v) {
        return Vector<Type>(st.ii() * v.x(), st.ii() * v.y(), st.ii() * v.z());
    }

    /**\brief Inner-product between a vector and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const Vector<Type> &v,
                                                          const SphericalTensor<Type> &st) {
        return Vector<Type>(v.x() * st.ii(), v.y() * st.ii(), v.z() * st.ii());
    }

    /**\brief Inner-product between a vector and a spherical tensor*/
    template <class Type>
    inline Vector<Type> operator&(const Vector<Type> &v, const SphericalTensor<Type> &st) {
        return Vector<Type>(v.x() * st.ii(), v.y() * st.ii(), v.z() * st.ii());
    }

    /**\brief Double-dot-product between a spherical tensor and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type operator&&(const SphericalTensor<Type> &st1,
                                                   const SphericalTensor<Type> &st2) noexcept {
        return 3 * st1.ii() * st2.ii();
    }

    /**\brief Division of a real by a sphericalTensor*/
    template <class Type>
    hur_nodiscard constexpr inline SphericalTensor<Type>
    operator/(const Type s, const SphericalTensor<Type> &st) {
        return SphericalTensor<Type>(s / st.ii());
    }

    template <class Type>
    hur_nodiscard constexpr inline Type magSqr(const SphericalTensor<Type> &st) noexcept {
        return 3 * magSqr(st.ii());
    }

    /**\brief Return the trace of a spherical tensor*/
    template <class Type> inline Type tr(const SphericalTensor<Type> &st) {
        return 3 * st.ii();
    }

    /**\brief Return the spherical part of a spherical tensor, i.e. itself*/
    template <class Type>
    hur_nodiscard constexpr inline SphericalTensor<Type>
    sph(const SphericalTensor<Type> &st) noexcept {
        return st;
    }

    /**\brief Return the determinant of a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type det(const SphericalTensor<Type> &st) noexcept {
        return st.ii() * st.ii() * st.ii();
    }

    /**\brief Return the inverse of a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SphericalTensor<Type> inv(const SphericalTensor<Type> &st) {
        return SphericalTensor<Type>(1.0 / st.ii());
    }

    /**
     * \brief The template class of identity tensor.
     */
    template <class Type> class Identity : public SphericalTensor<Type> {
    public:
        using Base = SphericalTensor<Type>;
        using value_type = typename Base::value_type;

    public:
        constexpr inline Identity() : Base(1) {}
        Identity(const Identity &other) = delete;
        Identity &operator=(const Identity &other) = delete;
        inline ~Identity() noexcept {}

        constexpr inline explicit operator integer() const noexcept { return 1; }
        constexpr inline explicit operator real() const noexcept { return 1.0; }
    };

    template <class Type> hur_nodiscard inline std::string toString(const Identity<Type> &) {
        std::string s;
        s.append("\n");
        s += "1, 0, 0";
        s.append("\n");
        s += "0, 1, 0";
        s.append("\n");
        s += "0, 0, 1";
        return s;
    }

} //  namespace OpenHurricane
