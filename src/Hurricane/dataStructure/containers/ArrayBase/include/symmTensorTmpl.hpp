/*!
 * \file symmTensorTmpl.hpp
 * \brief Header of Symm Tensor
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

#include "sphericalTensorTmpl.hpp"
#include "vectorSpace.hpp"
#include "tensorTmpl.hpp"

namespace OpenHurricane {
    /**
     * \brief The template class of SymmTensor.
     */
    template <class Type> class SymmTensor : public VectorSpace<Type, SymmTensor<Type>, 6> {
    public:
        using Base = VectorSpace<Type, SymmTensor<Type>, 6>;
        using value_type = typename Base::value_type;
        using tType = SymmTensor<value_type>;

        /**\brief Rank of SymmTensor is 2*/
        static constexpr int rank = 2;

        // Static data members

        static const SymmTensor I;

        /**\brief Component labeling enumeration*/
        enum components : short { XX, XY, XZ, YY, YZ, ZZ };

    public:
        /**\brief Construct null*/
        constexpr inline SymmTensor() : Base() {}

        /**\brief Construct initialized to zero*/
        constexpr inline SymmTensor(const OpenHurricane::zero) : Base(Zero) {}
        constexpr inline SymmTensor &operator=(const OpenHurricane::zero) {
            Base::operator=(Zero);
            return *this;
        }

        constexpr inline SymmTensor(const SymmTensor &other) : Base(other) {}
        constexpr inline SymmTensor &operator=(const SymmTensor &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }
        constexpr inline SymmTensor(const value_type &elem) : Base(elem) {}

        template <class otherType, class otherDerived>
        constexpr inline SymmTensor(const VectorSpace<otherType, otherDerived, 6> &other)
            : Base(other) {}

        constexpr inline SymmTensor(const SphericalTensor<value_type> &st) {
            this->operator[](XX) = st.ii();
            this->operator[](XY) = 0;
            this->operator[](XZ) = 0;
            this->operator[](YY) = st.ii();
            this->operator[](YZ) = 0;
            this->operator[](ZZ) = st.ii();
        }
        constexpr inline SymmTensor &operator=(const SphericalTensor<Type> &st) {
            this->operator[](XX) = st.ii();
            this->operator[](XY) = 0;
            this->operator[](XZ) = 0;
            this->operator[](YY) = st.ii();
            this->operator[](YZ) = 0;
            this->operator[](ZZ) = st.ii();
            return *this;
        }

        constexpr inline SymmTensor(const value_type txx, const value_type txy,
                                    const value_type txz, const value_type tyy,
                                    const value_type tyz, const value_type tzz) {
            this->operator[](XX) = txx;
            this->operator[](XY) = txy;
            this->operator[](XZ) = txz;
            this->operator[](YY) = tyy;
            this->operator[](YZ) = tyz;
            this->operator[](ZZ) = tzz;
        }

        inline ~SymmTensor() noexcept {}

        hur_nodiscard constexpr inline const value_type &xx() const noexcept {
            return this->operator[](XX);
        }
        hur_nodiscard constexpr inline const value_type &xy() const noexcept {
            return this->operator[](XY);
        }
        hur_nodiscard constexpr inline const value_type &xz() const noexcept {
            return this->operator[](XZ);
        }
        hur_nodiscard constexpr inline const value_type &yx() const noexcept {
            return this->operator[](XY);
        }
        hur_nodiscard constexpr inline const value_type &yy() const noexcept {
            return this->operator[](YY);
        }
        hur_nodiscard constexpr inline const value_type &yz() const noexcept {
            return this->operator[](YZ);
        }
        hur_nodiscard constexpr inline const value_type &zx() const noexcept {
            return this->operator[](XZ);
        }
        hur_nodiscard constexpr inline const value_type &zy() const noexcept {
            return this->operator[](YZ);
        }
        hur_nodiscard constexpr inline const value_type &zz() const noexcept {
            return this->operator[](ZZ);
        }

        hur_nodiscard constexpr inline value_type &xx() noexcept { return this->operator[](XX); }
        hur_nodiscard constexpr inline value_type &xy() noexcept { return this->operator[](XY); }
        hur_nodiscard constexpr inline value_type &xz() noexcept { return this->operator[](XZ); }
        hur_nodiscard constexpr inline value_type &yx() noexcept { return this->operator[](XY); }
        hur_nodiscard constexpr inline value_type &yy() noexcept { return this->operator[](YY); }
        hur_nodiscard constexpr inline value_type &yz() noexcept { return this->operator[](YZ); }
        hur_nodiscard constexpr inline value_type &zx() noexcept { return this->operator[](XZ); }
        hur_nodiscard constexpr inline value_type &zy() noexcept { return this->operator[](YZ); }
        hur_nodiscard constexpr inline value_type &zz() noexcept { return this->operator[](ZZ); }

        /**\brief Transpose*/
        hur_nodiscard constexpr inline const SymmTensor &transpose() const noexcept {
            return *this;
        }

        hur_nodiscard inline value_type magnitude() const {
            value_type Tsqr = sqr(xx()) + sqr(xy()) + sqr(xz()) + sqr(yx()) + sqr(yy()) +
                              sqr(yz()) + sqr(zx()) + sqr(zy()) + sqr(zz());
            return sqrt(Tsqr);
        }

        hur_nodiscard constexpr inline value_type magSqr() const noexcept {
            return sqr(xx()) + sqr(xy()) + sqr(xz()) + sqr(yx()) + sqr(yy()) + sqr(yz()) +
                   sqr(zx()) + sqr(zy()) + sqr(zz());
        }
    };

    /**\brief Product between a realand a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator*(const Type s,
                                                              const SymmTensor<Type> &st) {
        return SymmTensor<Type>(s * st.xx(), s * st.xy(), s * st.xz(), s * st.yy(), s * st.yz(),
                                s * st.zz());
    }

    /**\brief Hodge Dual operator (tensor -> vector)*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const SymmTensor<Type> &st) {
        return Vector<Type>(st.yz(), -st.xz(), st.xy());
    }

    /**\brief Inner-product between two symmetric tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator*(const SymmTensor<Type> &st1,
                                                          const SymmTensor<Type> &st2) {
        return Tensor<Type>(st1.xx() * st2.xx() + st1.xy() * st2.xy() + st1.xz() * st2.xz(),
                            st1.xx() * st2.xy() + st1.xy() * st2.yy() + st1.xz() * st2.yz(),
                            st1.xx() * st2.xz() + st1.xy() * st2.yz() + st1.xz() * st2.zz(),

                            st1.xy() * st2.xx() + st1.yy() * st2.xy() + st1.yz() * st2.xz(),
                            st1.xy() * st2.xy() + st1.yy() * st2.yy() + st1.yz() * st2.yz(),
                            st1.xy() * st2.xz() + st1.yy() * st2.yz() + st1.yz() * st2.zz(),

                            st1.xz() * st2.xx() + st1.yz() * st2.xy() + st1.zz() * st2.xz(),
                            st1.xz() * st2.xy() + st1.yz() * st2.yy() + st1.zz() * st2.yz(),
                            st1.xz() * st2.xz() + st1.yz() * st2.yz() + st1.zz() * st2.zz());
    }

    /**\brief Inner-product between two symmetric tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator&(const SymmTensor<Type> &st1,
                                                          const SymmTensor<Type> &st2) {
        return Tensor<Type>(st1.xx() * st2.xx() + st1.xy() * st2.xy() + st1.xz() * st2.xz(),
                            st1.xx() * st2.xy() + st1.xy() * st2.yy() + st1.xz() * st2.yz(),
                            st1.xx() * st2.xz() + st1.xy() * st2.yz() + st1.xz() * st2.zz(),

                            st1.xy() * st2.xx() + st1.yy() * st2.xy() + st1.yz() * st2.xz(),
                            st1.xy() * st2.xy() + st1.yy() * st2.yy() + st1.yz() * st2.yz(),
                            st1.xy() * st2.xz() + st1.yy() * st2.yz() + st1.yz() * st2.zz(),

                            st1.xz() * st2.xx() + st1.yz() * st2.xy() + st1.zz() * st2.xz(),
                            st1.xz() * st2.xy() + st1.yz() * st2.yy() + st1.zz() * st2.yz(),
                            st1.xz() * st2.xz() + st1.yz() * st2.yz() + st1.zz() * st2.zz());
    }

    /**\brief Double-dot-product between a symmetric tensor and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type operator&&(const SymmTensor<Type> &st1,
                                                   const SymmTensor<Type> &st2) {
        return (st1.xx() * st2.xx() + 2 * st1.xy() * st2.xy() + 2 * st1.xz() * st2.xz() +
                st1.yy() * st2.yy() + 2 * st1.yz() * st2.yz() + st1.zz() * st2.zz());
    }

    /**\brief Inner-product between a symmetric tensor and a vector*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const SymmTensor<Type> &st,
                                                          const Vector<Type> &v) {
        return Vector<Type>(st.xx() * v.x() + st.xy() * v.y() + st.xz() * v.z(),
                            st.xy() * v.x() + st.yy() * v.y() + st.yz() * v.z(),
                            st.xz() * v.x() + st.yz() * v.y() + st.zz() * v.z());
    }

    /**\brief Inner-product between a symmetric tensor and a vector*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator&(const SymmTensor<Type> &st,
                                                          const Vector<Type> &v) {
        return Vector<Type>(st.xx() * v.x() + st.xy() * v.y() + st.xz() * v.z(),
                            st.xy() * v.x() + st.yy() * v.y() + st.yz() * v.z(),
                            st.xz() * v.x() + st.yz() * v.y() + st.zz() * v.z());
    }

    /**\brief Inner-product between a vector and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const Vector<Type> &v,
                                                          const SymmTensor<Type> &st) {
        return Vector<Type>(v.x() * st.xx() + v.y() * st.xy() + v.z() * st.xz(),
                            v.x() * st.xy() + v.y() * st.yy() + v.z() * st.yz(),
                            v.x() * st.xz() + v.y() * st.yz() + v.z() * st.zz());
    }

    /**\brief Inner-product between a vector and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator&(const Vector<Type> &v,
                                                          const SymmTensor<Type> &st) {
        return Vector<Type>(v.x() * st.xx() + v.y() * st.xy() + v.z() * st.xz(),
                            v.x() * st.xy() + v.y() * st.yy() + v.z() * st.yz(),
                            v.x() * st.xz() + v.y() * st.yz() + v.z() * st.zz());
    }

    /**\brief Inner-sqr of a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> innerSqr(const SymmTensor<Type> &st) {
        return SymmTensor<Type>(st.xx() * st.xx() + st.xy() * st.xy() + st.xz() * st.xz(),
                                st.xx() * st.xy() + st.xy() * st.yy() + st.xz() * st.yz(),
                                st.xx() * st.xz() + st.xy() * st.yz() + st.xz() * st.zz(),

                                st.xy() * st.xy() + st.yy() * st.yy() + st.yz() * st.yz(),
                                st.xy() * st.xz() + st.yy() * st.yz() + st.yz() * st.zz(),

                                st.xz() * st.xz() + st.yz() * st.yz() + st.zz() * st.zz());
    }

    template <class Type> hur_nodiscard constexpr inline Type magSqr(const SymmTensor<Type> &st) {
        return (magSqr(st.xx()) + 2 * magSqr(st.xy()) + 2 * magSqr(st.xz()) + magSqr(st.yy()) +
                2 * magSqr(st.yz()) + magSqr(st.zz()));
    }

    /**\brief Return the trace of a symmetric tensor*/
    template <class Type> hur_nodiscard constexpr inline Type tr(const SymmTensor<Type> &st) {
        return st.xx() + st.yy() + st.zz();
    }

    /**\brief Return the spherical part of a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SphericalTensor<Type> sph(const SymmTensor<Type> &st) {
        return (1.0 / 3.0) * tr(st);
    }

    /**\brief Return the symmetric part of a symmetric tensor, i.e.itself*/
    template <class Type>
    hur_nodiscard constexpr inline const SymmTensor<Type> &symm(const SymmTensor<Type> &st) {
        return st;
    }

    /**\brief Return twice the symmetric part of a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> twoSymm(const SymmTensor<Type> &st) {
        return 2 * st;
    }

    /**\brief Return the deviatoric part of a symmetric tensor*/
    template <class Type> hur_nodiscard inline SymmTensor<Type> dev(const SymmTensor<Type> &st) {
        return st - SphericalTensor<Type>::oneThirdI * tr(st);
    }

    /**\brief Return the deviatoric part of a symmetric tensor*/
    template <class Type> hur_nodiscard inline SymmTensor<Type> dev2(const SymmTensor<Type> &st) {
        return st - SphericalTensor<Type>::twoThirdsI * tr(st);
    }

    /**\brief Return the determinant of a symmetric tensor*/
    template <class Type> hur_nodiscard constexpr inline Type det(const SymmTensor<Type> &st) {
        return (st.xx() * st.yy() * st.zz() + st.xy() * st.yz() * st.xz() +
                st.xz() * st.xy() * st.yz() - st.xx() * st.yz() * st.yz() -
                st.xy() * st.xy() * st.zz() - st.xz() * st.yy() * st.xz());
    }

    /**\brief Return the cofactor symmetric tensor of a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> cof(const SymmTensor<Type> &st) {
        return SymmTensor<Type>(
            st.yy() * st.zz() - st.yz() * st.yz(), st.xz() * st.yz() - st.xy() * st.zz(),
            st.xy() * st.yz() - st.xz() * st.yy(), st.xx() * st.zz() - st.xz() * st.xz(),
            st.xy() * st.xz() - st.xx() * st.yz(), st.xx() * st.yy() - st.xy() * st.xy());
    }

    /**\brief Return the inverse of a symmetric tensor give the determinant*/
    template <class Type>
    inline SymmTensor<Type> hur_nodiscard constexpr inv(const SymmTensor<Type> &st,
                                                        const Type detst) {
        return SymmTensor<Type>(
                   st.yy() * st.zz() - st.yz() * st.yz(), st.xz() * st.yz() - st.xy() * st.zz(),
                   st.xy() * st.yz() - st.xz() * st.yy(), st.xx() * st.zz() - st.xz() * st.xz(),
                   st.xy() * st.xz() - st.xx() * st.yz(), st.xx() * st.yy() - st.xy() * st.xy()) /
               detst;
    }

    /**\brief Return the inverse of a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> inv(const SymmTensor<Type> &st) {
        return inv(st, det(st));
    }

    /**\brief Return the 1st invariant of a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type invariantI(const SymmTensor<Type> &st) {
        return tr(st);
    }

    /**\brief Return the 2nd invariant of a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type invariantII(const SymmTensor<Type> &st) {
        return (st.xx() * st.yy() + st.yy() * st.zz() + st.xx() * st.zz() - sqr(st.xy()) -
                sqr(st.yz()) - sqr(st.xz()));
    }

    /**\brief Return the 3rd invariant of a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type invariantIII(const SymmTensor<Type> &st) {
        return det(st);
    }

    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator+(const SphericalTensor<Type> &spt1,
                                                              const SymmTensor<Type> &st2) {
        return SymmTensor<Type>(spt1.ii() + st2.xx(), st2.xy(), st2.xz(), spt1.ii() + st2.yy(),
                                st2.yz(), spt1.ii() + st2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator+(const SymmTensor<Type> &st1,
                                                              const SphericalTensor<Type> &spt2) {
        return SymmTensor<Type>(st1.xx() + spt2.ii(), st1.xy(), st1.xz(), st1.yy() + spt2.ii(),
                                st1.yz(), st1.zz() + spt2.ii());
    }

    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator-(const SphericalTensor<Type> &spt1,
                                                              const SymmTensor<Type> &st2) {
        return SymmTensor<Type>(spt1.ii() - st2.xx(), -st2.xy(), -st2.xz(), spt1.ii() - st2.yy(),
                                -st2.yz(), spt1.ii() - st2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator-(const SymmTensor<Type> &st1,
                                                              const SphericalTensor<Type> &spt2) {
        return SymmTensor<Type>(st1.xx() - spt2.ii(), st1.xy(), st1.xz(), st1.yy() - spt2.ii(),
                                st1.yz(), st1.zz() - spt2.ii());
    }

    /**\brief Inner-product between a spherical symmetric tensor and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator*(const SphericalTensor<Type> &spt1,
                                                              const SymmTensor<Type> &st2) {
        return SymmTensor<Type>(spt1.ii() * st2.xx(), spt1.ii() * st2.xy(), spt1.ii() * st2.xz(),
                                spt1.ii() * st2.yy(), spt1.ii() * st2.yz(), spt1.ii() * st2.zz());
    }

    /**\brief Inner-product between a spherical symmetric tensor and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator&(const SphericalTensor<Type> &spt1,
                                                              const SymmTensor<Type> &st2) {
        return SymmTensor<Type>(spt1.ii() * st2.xx(), spt1.ii() * st2.xy(), spt1.ii() * st2.xz(),
                                spt1.ii() * st2.yy(), spt1.ii() * st2.yz(), spt1.ii() * st2.zz());
    }

    /**\brief Inner-product between a symmetric tensor and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator*(const SymmTensor<Type> &st1,
                                                              const SphericalTensor<Type> &spt2) {
        return SymmTensor<Type>(st1.xx() * spt2.ii(), st1.xy() * spt2.ii(), st1.xz() * spt2.ii(),
                                st1.yy() * spt2.ii(), st1.yz() * spt2.ii(), st1.zz() * spt2.ii());
    }

    /**\brief Inner-product between a symmetric tensor and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> operator&(const SymmTensor<Type> &st1,
                                                              const SphericalTensor<Type> &spt2) {
        return SymmTensor<Type>(st1.xx() * spt2.ii(), st1.xy() * spt2.ii(), st1.xz() * spt2.ii(),
                                st1.yy() * spt2.ii(), st1.yz() * spt2.ii(), st1.zz() * spt2.ii());
    }

    /**\brief Double-dot-product between a spherical tensor and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type operator&&(const SphericalTensor<Type> &spt1,
                                                   const SymmTensor<Type> &st2) {
        return (spt1.ii() * st2.xx() + spt1.ii() * st2.yy() + spt1.ii() * st2.zz());
    }

    /**\brief Double-dot-product between a symmetric tensor and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type operator&&(const SymmTensor<Type> &st1,
                                                   const SphericalTensor<Type> &spt2) {
        return (st1.xx() * spt2.ii() + st1.yy() * spt2.ii() + st1.zz() * spt2.ii());
    }

    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> sqr(const Vector<Type> &v) {
        return SymmTensor<Type>(v.x() * v.x(), v.x() * v.y(), v.x() * v.z(), v.y() * v.y(),
                                v.y() * v.z(), v.z() * v.z());
    }

} //  namespace OpenHurricane
