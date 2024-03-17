/*!
 * \file diagTensorTmpl.hpp
 * \brief Header of diagonal Tensor
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

#include "sphericalTensorTmpl.hpp"
#include "symmTensorTmpl.hpp"
#include "tensorTmpl.hpp"

namespace OpenHurricane {

    /**
     * \brief The template class of DiagTensor.
     */
    template <class Type> class DiagTensor : public VectorSpace<Type, DiagTensor<Type>, 3> {
    public:
        using Base = VectorSpace<Type, DiagTensor<Type>, 3>;
        using value_type = typename Base::value_type;
        using tType = DiagTensor<value_type>;

        /*!\brief Rank of DiagTensor is 2*/
        static constexpr int rank = 2;

        static const DiagTensor I;

        /*!\brief Component labeling enumeration*/
        enum components : short { XX, YY, ZZ };

    public:
        /*!\brief Construct null*/
        constexpr inline DiagTensor() : Base() {}
        constexpr inline DiagTensor(const OpenHurricane::zero) : Base(Zero) {}
        constexpr inline DiagTensor &operator=(const OpenHurricane::zero) {
            Base::operator=(Zero);
            return *this;
        }

        constexpr inline DiagTensor(const DiagTensor &other) : Base(other) {}
        constexpr inline DiagTensor &operator=(const DiagTensor &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }
        constexpr inline DiagTensor(const value_type &elem) : Base(elem) {}

        template <class otherType, class otherDerived>
        constexpr inline DiagTensor(const VectorSpace<otherType, otherDerived, 3> &other)
            : Base(other) {}

        /*!\brief Construct given three components*/
        constexpr inline DiagTensor(const value_type &txx, const value_type &tyy,
                                    const value_type &tzz) {
            this->operator[](XX) = txx;
            this->operator[](YY) = tyy;
            this->operator[](ZZ) = tzz;
        }

        hur_nodiscard constexpr inline const value_type &xx() const noexcept {
            return this->operator[](XX);
        }
        hur_nodiscard constexpr inline const value_type &yy() const noexcept {
            return this->operator[](YY);
        }
        hur_nodiscard constexpr inline const value_type &zz() const noexcept {
            return this->operator[](ZZ);
        }

        hur_nodiscard constexpr inline value_type &xx() noexcept { return this->operator[](XX); }
        hur_nodiscard constexpr inline value_type &yy() noexcept { return this->operator[](YY); }
        hur_nodiscard constexpr inline value_type &zz() noexcept { return this->operator[](ZZ); }

        /**\brief Transpose*/
        hur_nodiscard constexpr inline const DiagTensor &transpose() const noexcept {
            return *this;
        }

        hur_nodiscard inline value_type magnitude() const {
            value_type Tsqr = sqr(xx()) + sqr(yy()) + sqr(zz());
            return sqrt(Tsqr);
        }

        hur_nodiscard constexpr inline value_type magSqr() const noexcept {
            return sqr(xx()) + sqr(yy()) + sqr(zz());
        }
    };

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator+(const DiagTensor<Type> &dt1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(dt1.xx() + t2.xx(), t2.xy(), t2.xz(), t2.yx(), dt1.yy() + t2.yy(),
                            t2.yz(), t2.zx(), t2.zy(), dt1.zz() + t2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator+(const Tensor<Type> &t1,
                                                          const DiagTensor<Type> &dt2) {
        return Tensor<Type>(t1.xx() + dt2.xx(), t1.xy(), t1.xz(), t1.yx(), t1.yy() + dt2.yy(),
                            t1.yz(), t1.zx(), t1.zy(), t1.zz() + dt2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator+(const DiagTensor<Type> &dt1,
                                                              const SymmTensor<Type> &st2) {
        return SymmTensor<Type>(dt1.xx() + st2.xx(), st2.xy(), st2.xz(), dt1.yy() + st2.yy(),
                                st2.yz(), dt1.zz() + st2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator+(const SymmTensor<Type> &st1,
                                                              const DiagTensor<Type> &dt2) {
        return SymmTensor<Type>(st1.xx() + dt2.xx(), st1.xy(), st1.xz(), st1.yy() + dt2.yy(),
                                st1.yz(), st1.zz() + dt2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator-(const DiagTensor<Type> &dt1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(dt1.xx() - t2.xx(), -t2.xy(), -t2.xz(), -t2.yx(), dt1.yy() - t2.yy(),
                            -t2.yz(), -t2.zx(), -t2.zy(), dt1.zz() - t2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator-(const Tensor<Type> &t1,
                                                          const DiagTensor<Type> &dt2) {
        return Tensor<Type>(t1.xx() - dt2.xx(), t1.xy(), t1.xz(), t1.yx(), t1.yy() - dt2.yy(),
                            t1.yz(), t1.zx(), t1.zy(), t1.zz() - dt2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator-(const DiagTensor<Type> &dt1,
                                                              const SymmTensor<Type> &st2) {
        return SymmTensor<Type>(dt1.xx() - st2.xx(), -st2.xy(), -st2.xz(), dt1.yy() - st2.yy(),
                                -st2.yz(), dt1.zz() - st2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator-(const SymmTensor<Type> &st1,
                                                              const DiagTensor<Type> &dt2) {
        return SymmTensor<Type>(st1.xx() - dt2.xx(), st1.xy(), st1.xz(), st1.yy() - dt2.yy(),
                                st1.yz(), st1.zz() - dt2.zz());
    }

    /**\brief Inner-product between two diagonal tensors */
    template <class Type>
    hur_nodiscard constexpr inline DiagTensor<Type> operator&(const DiagTensor<Type> &dt1,
                                                              const DiagTensor<Type> &dt2) {
        return DiagTensor<Type>(dt1.xx() * dt2.xx(), dt1.yy() * dt2.yy(), dt1.zz() * dt2.zz());
    }

    /**\brief Inner-product between two diagonal tensors */
    template <class Type>
    hur_nodiscard constexpr inline DiagTensor<Type> operator*(const DiagTensor<Type> &dt1,
                                                              const DiagTensor<Type> &dt2) {
        return DiagTensor<Type>(dt1.xx() * dt2.xx(), dt1.yy() * dt2.yy(), dt1.zz() * dt2.zz());
    }

    /**\brief Inner-product between a diagonal tensor and a tensor */
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator&(const DiagTensor<Type> &dt1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(dt1.xx() * t2.xx(), dt1.xx() * t2.xy(), dt1.xx() * t2.xz(),
                            dt1.yy() * t2.yx(), dt1.yy() * t2.yy(), dt1.yy() * t2.yz(),
                            dt1.zz() * t2.zx(), dt1.zz() * t2.zy(), dt1.zz() * t2.zz());
    }

    /**\brief Inner-product between a diagonal tensor and a tensor */
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator*(const DiagTensor<Type> &dt1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(dt1.xx() * t2.xx(), dt1.xx() * t2.xy(), dt1.xx() * t2.xz(),
                            dt1.yy() * t2.yx(), dt1.yy() * t2.yy(), dt1.yy() * t2.yz(),
                            dt1.zz() * t2.zx(), dt1.zz() * t2.zy(), dt1.zz() * t2.zz());
    }

    /**\brief Inner-product between a tensor and a diagonal tensor */
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator&(const Tensor<Type> &t1,
                                                          const DiagTensor<Type> &dt2) {
        return Tensor<Type>(t1.xx() * dt2.xx(), t1.xy() * dt2.yy(), t1.xz() * dt2.zz(),
                            t1.yx() * dt2.xx(), t1.yy() * dt2.yy(), t1.yz() * dt2.zz(),
                            t1.zx() * dt2.xx(), t1.zy() * dt2.yy(), t1.zz() * dt2.zz());
    }

    /**\brief Inner-product between a tensor and a diagonal tensor */
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator*(const Tensor<Type> &t1,
                                                          const DiagTensor<Type> &dt2) {
        return Tensor<Type>(t1.xx() * dt2.xx(), t1.xy() * dt2.yy(), t1.xz() * dt2.zz(),
                            t1.yx() * dt2.xx(), t1.yy() * dt2.yy(), t1.yz() * dt2.zz(),
                            t1.zx() * dt2.xx(), t1.zy() * dt2.yy(), t1.zz() * dt2.zz());
    }

    /**\brief Inner-product between a diagonal tensor and a vector */
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator&(const DiagTensor<Type> &dt,
                                                          const Vector<Type> &v) {
        return Vector<Type>(dt.xx() * v.x(), dt.yy() * v.y(), dt.zz() * v.z());
    }

    /**\brief Inner-product between a diagonal tensor and a vector */
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const DiagTensor<Type> &dt,
                                                          const Vector<Type> &v) {
        return Vector<Type>(dt.xx() * v.x(), dt.yy() * v.y(), dt.zz() * v.z());
    }

    /**\brief Inner-product between a vector and a diagonal tensor */
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator&(const Vector<Type> &v,
                                                          const DiagTensor<Type> &dt) {
        return Vector<Type>(v.x() * dt.xx(), v.y() * dt.yy(), v.z() * dt.zz());
    }

    /**\brief Inner-product between a vector and a diagonal tensor */
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const Vector<Type> &v,
                                                          const DiagTensor<Type> &dt) {
        return Vector<Type>(v.x() * dt.xx(), v.y() * dt.yy(), v.z() * dt.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline DiagTensor<Type> operator*(const Type s,
                                                              const DiagTensor<Type> &dt) {
        return DiagTensor<Type>(s * dt.xx(), s * dt.yy(), s * dt.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline DiagTensor<Type> operator*(const DiagTensor<Type> &dt,
                                                              const Type s) {
        return DiagTensor<Type>(s * dt.xx(), s * dt.yy(), s * dt.zz());
    }

    /**\brief Division of a real by a diagonalTensor*/
    template <class Type>
    hur_nodiscard constexpr inline DiagTensor<Type> operator/(const Type s,
                                                              const DiagTensor<Type> &dt) {
        return DiagTensor<Type>(s / dt.xx(), s / dt.yy(), s / dt.zz());
    }

    /**\brief Division of a vector by a diagonalTensor*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator/(const Vector<Type> v,
                                                          const DiagTensor<Type> &dt) {
        return Vector<Type>(v.x() / dt.xx(), v.y() / dt.yy(), v.z() / dt.zz());
    }

    /**\brief Return the trace of a diagonal tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type tr(const DiagTensor<Type> &dt) noexcept {
        return dt.xx() + dt.yy() + dt.zz();
    }

    /**\brief Return the spherical part of a diagonal tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SphericalTensor<Type> sph(const DiagTensor<Type> &dt) noexcept {
        return 0.5 * tr(dt);
    }

    /**\brief Return the determinant of a diagonal tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type det(const DiagTensor<Type> &t) noexcept {
        return t.xx() * t.yy() * t.zz();
    }

    /**\brief Return the inverse of a diagonal tensor*/
    template <class Type>
    hur_nodiscard constexpr inline DiagTensor<Type> inv(const DiagTensor<Type> &dt) {
        return DiagTensor<Type>(1.0 / dt.xx(), 1.0 / dt.yy(), 1.0 / dt.zz());
    }

    /**\brief Return the diagonal of a symmetric tensor as a diagonal tensor*/
    template <class Type>
    hur_nodiscard constexpr inline DiagTensor<Type> diag(const SymmTensor<Type> &t) {
        return DiagTensor<Type>(t.xx(), t.yy(), t.zz());
    }

    /**\brief Return the diagonal of a tensor as a diagonal tensor*/
    template <class Type>
    hur_nodiscard constexpr inline DiagTensor<Type> diag(const Tensor<Type> &t) {
        return DiagTensor<Type>(t.xx(), t.yy(), t.zz());
    }
} //  namespace OpenHurricane
