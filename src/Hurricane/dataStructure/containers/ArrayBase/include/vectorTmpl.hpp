/*!
 * \file vectorTmpl.hpp
 * \brief Header of vector3D.
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
#include "vectorSpace.hpp"
namespace OpenHurricane {
    /**
     *\brief The template class of vector.
     *        - Rank of this vector = 1.
     *        - Dimensional sets of this vector = 3.
     */
    template <class Type> class Vector : public VectorSpace<Type, Vector<Type>, 3> {
    public:
        using Base = VectorSpace<Type, Vector<Type>, 3>;
        using value_type = typename Base::value_type;
        using vType = Vector<Type>;

        static constexpr int rank = 1;

        /**\brief Component labeling enumeration.*/
        enum components : short { X = 0, Y = 1, Z = 2 };

    public:
        constexpr inline Vector() : Base() {}

        constexpr inline Vector(const OpenHurricane::zero) : Base(Zero) {}
        constexpr inline Vector &operator=(const OpenHurricane::zero) {
            Base::operator=(Zero);
            return *this;
        }

        constexpr inline Vector(const Vector &other) : Base(other) {}
        constexpr inline Vector &operator=(const Vector &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        constexpr inline Vector(const value_type &elem) : Base(elem) {}

        template <class otherType, class otherDerived>
        constexpr inline Vector(const VectorSpace<otherType, otherDerived, 3> &other)
            : Base(other) {}

        constexpr inline Vector(const value_type &vx, const value_type &vy, const value_type &vz) {
            this->operator[](X) = vx;
            this->operator[](Y) = vy;
            this->operator[](Z) = vz;
        }

        inline Vector(IStringStream &is) : Base(is) {}
        inline ~Vector() noexcept {}

        hur_nodiscard constexpr inline const value_type &x() const noexcept {
            return this->operator[](X);
        }
        hur_nodiscard constexpr inline const value_type &y() const noexcept {
            return this->operator[](Y);
        }
        hur_nodiscard constexpr inline const value_type &z() const noexcept {
            return this->operator[](Z);
        }

        hur_nodiscard constexpr inline value_type &x() noexcept { return this->operator[](X); }
        hur_nodiscard constexpr inline value_type &y() noexcept { return this->operator[](Y); }
        hur_nodiscard constexpr inline value_type &z() noexcept { return this->operator[](Z); }

        /*!\brief Return a normalized vector of Vector*/
        hur_nodiscard inline Vector normalized() const {
#ifdef HUR_DEBUG
            if (this->operator==(OpenHurricane::Zero)) {
                LFatal("Attempt to normalize a zero vector!");
            }
#endif // HUR_DEBUG
            value_type mag = Base::magnitude();
            vType vt(this->operator[](X) / mag, this->operator[](Y) / mag,
                     this->operator[](Z) / mag);
            return vt;
        }

        /*!\brief Vector projection
         *  Return the projection vector on Vector n
         */
        hur_nodiscard inline vType vProject(const Vector &n) const {
            return ((*this) * n.normalized()) * n.normalized();
        }

        /*!\brief Real projection
         *  Return the projection real on Vector n
         */
        hur_nodiscard inline value_type sProject(const Vector &n) const {
            return ((*this) * n.normalized());
        }

        constexpr inline void operator+=(const Vector &v) noexcept { Base::operator+=(v); }
        constexpr inline void operator-=(const Vector &v) noexcept { Base::operator-=(v); }

        constexpr inline void operator*=(const value_type &s) noexcept { Base::operator*=(s); }

        constexpr inline void operator*=(const Vector &v) noexcept {
            for (int i = 0; i < 3; ++i) {
                this->operator[](i) *= v.v_[i];
            }
        }
        constexpr inline void operator/=(const value_type &s) { Base::operator/=(s); }

        hur_nodiscard constexpr inline bool operator==(const Vector &other) const noexcept {
            return Base::operator==(other);
        }
        hur_nodiscard constexpr inline bool operator!=(const Vector &other) const noexcept {
            return !this->operator==(other);
        }

        hur_nodiscard constexpr inline bool operator==(const OpenHurricane::zero) const noexcept {
            return Base::operator==(Zero);
        }
        hur_nodiscard constexpr inline bool operator!=(const OpenHurricane::zero) const noexcept {
            return !this->operator==(Zero);
        }
    };

    template <class Type> class rankType<Type, 1> {
    public:
        using type = Vector<Type>;
    };

    template <class Type> class innerProduct<Vector<Type>, real> {
    public:
        using type = real;
    };

    /*!\brief Dot product of two vectors. (inner-products)
     *  Return a real
     *  (x1,y1,z1)*(x2,y2,z2)^T=x1x2+y1y2+z1z2
     */
    template <class Type>
    hur_nodiscard constexpr hur_forceinline Type operator*(const Vector<Type> &v1,
                                                           const Vector<Type> &v2) {
        return (v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z());
    }

    /*!\brief Cross-product of two vectors
     *  Return a vector
     *  (x1,y1,z1)^(x2,y2,z2) = (y1z2-z1y2, z1x2-x1z2, x1y2-y1x2)
     */
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator^(const Vector<Type> &v1,
                                                          const Vector<Type> &v2) {
        return Vector<Type>((v1.y() * v2.z() - v1.z() * v2.y()),
                            (v1.z() * v2.x() - v1.x() * v2.z()),
                            (v1.x() * v2.y() - v1.y() * v2.x()));
    }

    /*!\brief Real multiplication of vectors
     *  Return a vector
     */
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const Type &k, const Vector<Type> &v) {
        return Vector<Type>(k * v.x(), k * v.y(), k * v.z());
    }

    /*!\brief Real multiplication of vectors
     *  Return a vector
     */
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const Vector<Type> &v, const Type &k) {
        return Vector<Type>(k * v.x(), k * v.y(), k * v.z());
    }

    /*!\brief Real multiplication of vectors (outer product)
     *  Return a vector
     */
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator&(const Type &k, const Vector<Type> &v) {
        return Vector<Type>(k * v.x(), k * v.y(), k * v.z());
    }

    /*!\brief Real multiplication of vectors (outer product)
     *  Return a vector
     */
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator&(const Vector<Type> &v, const Type &k) {
        Vector<Type> tmp;
        for (int i = 0; i < 3; ++i) {
            tmp[i] = k * v[i];
        }
        return tmp;
    }

    /*!\brief Return a vector*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator/(const Vector<Type> &v, const Type &k) {
#ifdef HUR_DEBUG
        if (k == 0) {
            LFatal("divided by zero!");
        }
#endif // HUR_DEBUG
        return Vector<Type>(v.x() / k, v.y() / k, v.z() / k);
    }

    /*!
     * \brief Negating a vector
     * \retval A vector
     */
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator-(const Vector<Type> &v) {
        return Vector<Type>(-v.x(), -v.y(), -v.z());
    }

    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator-(const Vector<Type> &v1,
                                                          const Vector<Type> &v2) {
        return Vector<Type>((v1.x() - v2.x()), (v1.y() - v2.y()), (v1.z() - v2.z()));
    }

    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator+(const Vector<Type> &v1,
                                                          const Vector<Type> &v2) {
        return Vector<Type>((v1.x() + v2.x()), (v1.y() + v2.y()), (v1.z() + v2.z()));
    }

    /*!\brief included angle of vectors
     *  Return a real
     */
    template <class Type>
    hur_nodiscard inline Type cos(const Vector<Type> &v1, const Vector<Type> &v2) {
#ifdef HUR_DEBUG
        if ((v1 == OpenHurricane::Zero) && (v2 == OpenHurricane::Zero)) {
            LFatal("attempt to calculate the included angle of two zero vectors!");
        }
#endif // HUR_DEBUG

        return ((v1 * v2) / (v1.magnitude() * v2.magnitude()));
    }

    /*!\brief Distance between two points (vectors)
     *  Return a real
     */
    template <class Type>
    hur_nodiscard inline Type dist(const Vector<Type> &v1, const Vector<Type> &v2) {
        Type dx = v1.x() - v2.x();
        Type dy = v1.y() - v2.y();
        Type dz = v1.z() - v2.z();
        return (sqrt(dx * dx + dy * dy + dz * dz));
    }

    /*!\brief The magnitude of Vector
     *  Return a real
     */
    template <class Type> hur_nodiscard inline Type mag(const Vector<Type> &v1) {
        return (sqrt(v1.x() * v1.x() + v1.y() * v1.y() + v1.z() * v1.z()));
    }

    /*!\brief The magnitude square of Vector
     *  Return a real
     */
    template <class Type>
    hur_nodiscard constexpr inline Type magSqr(const Vector<Type> &v1) noexcept {
        return (v1.x() * v1.x() + v1.y() * v1.y() + v1.z() * v1.z());
    }

    template <class Type> hur_nodiscard constexpr inline Type div(const Vector<Type> &v1) noexcept {
        return v1.x() + v1.y() + v1.z();
    }
} // namespace OpenHurricane
