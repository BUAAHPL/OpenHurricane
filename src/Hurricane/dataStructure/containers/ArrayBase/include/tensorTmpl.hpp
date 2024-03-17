/*!
 * \file tensorTmpl.hpp
 * \brief Header of Tensor
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
#include "vectorTmpl.hpp"

namespace OpenHurricane {

    template <class Type> class SymmTensor;
    template <class Type> class DiagTensor;

    /**
     *\brief The template class of tensor.
     *        - Rank of this tensor = 2.
     */
    template <class Type> class Tensor : public VectorSpace<Type, Tensor<Type>, 9> {
    public:
        using Base = VectorSpace<Type, Tensor<Type>, 9>;
        using value_type = typename Base::value_type;
        using tType = Tensor<value_type>;

        /**\brief Rank of Tensor is 2*/
        static constexpr int rank = 2;

        static const Tensor I;

        /**\brief Component labeling enumeration*/
        enum components : short { XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };

    public:
        /**\brief Construct null*/
        constexpr inline Tensor() : Base() {}

        constexpr inline Tensor(const OpenHurricane::zero) : Base(Zero) {}
        constexpr inline Tensor &operator=(const OpenHurricane::zero) {
            Base::operator=(Zero);
            return *this;
        }

        constexpr inline Tensor(const Tensor &other) : Base(other) {}
        constexpr inline Tensor &operator=(const Tensor &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }
        constexpr inline Tensor(const value_type &elem) : Base(elem) {}

        template <class otherType, class otherDerived>
        constexpr inline Tensor(const VectorSpace<otherType, otherDerived, 9> &other)
            : Base(other) {}

        /**\brief Construct given SphericalTensor*/
        constexpr inline Tensor(const SphericalTensor<value_type> &st) {
            this->operator[](XX) = st.ii();
            this->operator[](XY) = 0;
            this->operator[](XZ) = 0;
            this->operator[](YX) = 0;
            this->operator[](YY) = st.ii();
            this->operator[](YZ) = 0;
            this->operator[](ZX) = 0;
            this->operator[](ZY) = 0;
            this->operator[](ZZ) = st.ii();
        }

        /**\brief Construct given SymmTensor*/
        constexpr inline Tensor(const SymmTensor<value_type> &st) {
            this->operator[](XX) = st.xx();
            this->operator[](XY) = st.xy();
            this->operator[](XZ) = st.xz();
            this->operator[](YX) = st.xy();
            this->operator[](YY) = st.yy();
            this->operator[](YZ) = st.yz();
            this->operator[](ZX) = st.xz();
            this->operator[](ZY) = st.yz();
            this->operator[](ZZ) = st.zz();
        }

        /**\brief Construct given triad*/
        constexpr inline Tensor(const Vector<Vector<value_type>> &tr) {
            this->operator[](XX) = tr.x().x();
            this->operator[](XY) = tr.x().y();
            this->operator[](XZ) = tr.x().z();

            this->operator[](YX) = tr.y().x();
            this->operator[](YY) = tr.y().y();
            this->operator[](YZ) = tr.y().z();

            this->operator[](ZX) = tr.z().x();
            this->operator[](ZY) = tr.z().y();
            this->operator[](ZZ) = tr.z().z();
        }

        /**\brief Construct given the three vector components*/
        constexpr inline Tensor(const Vector<value_type> &x, const Vector<value_type> &y,
                                const Vector<value_type> &z) {
            this->operator[](XX) = x.x();
            this->operator[](XY) = x.y();
            this->operator[](XZ) = x.z();
            this->operator[](YX) = y.x();
            this->operator[](YY) = y.y();
            this->operator[](YZ) = y.z();
            this->operator[](ZX) = z.x();
            this->operator[](ZY) = z.y();
            this->operator[](ZZ) = z.z();
        }

        /**\brief Construct given the nine components*/
        constexpr inline Tensor(const value_type txx, const value_type txy, const value_type txz,
                                const value_type tyx, const value_type tyy, const value_type tyz,
                                const value_type tzx, const value_type tzy, const value_type tzz) {
            this->operator[](XX) = txx;
            this->operator[](XY) = txy;
            this->operator[](XZ) = txz;
            this->operator[](YX) = tyx;
            this->operator[](YY) = tyy;
            this->operator[](YZ) = tyz;
            this->operator[](ZX) = tzx;
            this->operator[](ZY) = tzy;
            this->operator[](ZZ) = tzz;
        }

        inline ~Tensor() noexcept {}

        hur_nodiscard constexpr inline Vector<value_type> x() const {
            return Vector<value_type>(this->operator[](XX), this->operator[](XY),
                                      this->operator[](XZ));
        }

        hur_nodiscard constexpr inline Vector<value_type> y() const {
            return Vector<value_type>(this->operator[](YX), this->operator[](YY),
                                      this->operator[](YZ));
        }

        hur_nodiscard constexpr inline Vector<value_type> z() const {
            return Vector<value_type>(this->operator[](ZX), this->operator[](ZY),
                                      this->operator[](ZZ));
        }

        hur_nodiscard constexpr inline Vector<value_type> d() const {
            return Vector<value_type>(this->operator[](XX), this->operator[](YY),
                                      this->operator[](ZZ));
        }

        hur_nodiscard inline Vector<value_type> vectorComponent(const int id) const {
            switch (id) {
            case 0:
                return x();
                break;
            case 1:
                return y();
                break;
            case 2:
                return z();
                break;
            default:
                LFatal("invalid index: %d", id);
                return x();
                break;
            }
        }

        // Const element access functions for a 3x3
        // Compile-time errors are generated for inappropriate use

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
            return this->operator[](YX);
        }
        hur_nodiscard constexpr inline const value_type &yy() const noexcept {
            return this->operator[](YY);
        }
        hur_nodiscard constexpr inline const value_type &yz() const noexcept {
            return this->operator[](YZ);
        }
        hur_nodiscard constexpr inline const value_type &zx() const noexcept {
            return this->operator[](ZX);
        }
        hur_nodiscard constexpr inline const value_type &zy() const noexcept {
            return this->operator[](ZY);
        }
        hur_nodiscard constexpr inline const value_type &zz() const noexcept {
            return this->operator[](ZZ);
        }
        hur_nodiscard constexpr inline value_type &xx() noexcept { return this->operator[](XX); }
        hur_nodiscard constexpr inline value_type &xy() noexcept { return this->operator[](XY); }
        hur_nodiscard constexpr inline value_type &xz() noexcept { return this->operator[](XZ); }
        hur_nodiscard constexpr inline value_type &yx() noexcept { return this->operator[](YX); }
        hur_nodiscard constexpr inline value_type &yy() noexcept { return this->operator[](YY); }
        hur_nodiscard constexpr inline value_type &yz() noexcept { return this->operator[](YZ); }
        hur_nodiscard constexpr inline value_type &zx() noexcept { return this->operator[](ZX); }
        hur_nodiscard constexpr inline value_type &zy() noexcept { return this->operator[](ZY); }
        hur_nodiscard constexpr inline value_type &zz() noexcept { return this->operator[](ZZ); }

        hur_nodiscard constexpr inline const value_type &operator()(const int &i,
                                                                    const int &j) const noexcept {
#ifdef HUR_DEBUG
            if (i > 2 || j > 2) {
                LFatal("indices out of range");
            }
#endif
            return this->operator[](i * 3 + j);
        }

        // (i, j) element access operator
        hur_nodiscard constexpr inline value_type &operator()(const int &i, const int &j) noexcept {
#ifdef HUR_DEBUG
            if (i > 2 || j > 2) {
                LFatal("indices out of range");
            }
#endif
            return this->operator[](i * 3 + j);
        }

        /**\brief Return transpose*/
        hur_nodiscard constexpr inline Tensor transpose() const {
            return tType(xx(), yx(), zx(), xy(), yy(), zy(), xz(), yz(), zz());
        }

        /**\brief Return inverse*/
        hur_nodiscard constexpr inline Tensor inv() const { return OpenHurricane::inv(*this); }

        /**\brief Inner-product with a Tensor*/
        constexpr inline void operator*=(const Tensor<Type> &t) {
            *this = (Tensor<Type>(this->xx() * t.xx() + this->xy() * t.yx() + this->xz() * t.zx(),
                                  this->xx() * t.xy() + this->xy() * t.yy() + this->xz() * t.zy(),
                                  this->xx() * t.xz() + this->xy() * t.yz() + this->xz() * t.zz(),

                                  this->yx() * t.xx() + this->yy() * t.yx() + this->yz() * t.zx(),
                                  this->yx() * t.xy() + this->yy() * t.yy() + this->yz() * t.zy(),
                                  this->yx() * t.xz() + this->yy() * t.yz() + this->yz() * t.zz(),

                                  this->zx() * t.xx() + this->zy() * t.yx() + this->zz() * t.zx(),
                                  this->zx() * t.xy() + this->zy() * t.yy() + this->zz() * t.zy(),
                                  this->zx() * t.xz() + this->zy() * t.yz() + this->zz() * t.zz()));
        }

        /**\brief Assign to an equivalent vector space*/
        template <class otherType, class otherDerived>
        constexpr inline void operator=(const VectorSpace<otherType, otherDerived, 9> &vs) {
            Base::operator=(vs);
        }

        constexpr inline Tensor &operator=(const SphericalTensor<Type> &st) {
            this->operator[](XX) = st.ii();
            this->operator[](XY) = 0;
            this->operator[](XZ) = 0;
            this->operator[](YX) = 0;
            this->operator[](YY) = st.ii();
            this->operator[](YZ) = 0;
            this->operator[](ZX) = 0;
            this->operator[](ZY) = 0;
            this->operator[](ZZ) = st.ii();
            return *this;
        }

        constexpr inline Tensor &operator=(const SymmTensor<Type> &st) {
            this->operator[](XX) = st.xx();
            this->operator[](XY) = st.xy();
            this->operator[](XZ) = st.xz();
            this->operator[](YX) = st.xy();
            this->operator[](YY) = st.yy();
            this->operator[](YZ) = st.yz();
            this->operator[](ZX) = st.xz();
            this->operator[](ZY) = st.yz();
            this->operator[](ZZ) = st.zz();
            return *this;
        }

        constexpr inline Tensor &operator=(const Vector<Vector<Type>> &tr) {
            this->operator[](XX) = tr.x().x();
            this->operator[](XY) = tr.x().y();
            this->operator[](XZ) = tr.x().z();

            this->operator[](YX) = tr.y().x();
            this->operator[](YY) = tr.y().y();
            this->operator[](YZ) = tr.y().z();

            this->operator[](ZX) = tr.z().x();
            this->operator[](ZY) = tr.z().y();
            this->operator[](ZZ) = tr.z().z();
            return *this;
        }
    };

    template <class Type> class rankType<Type, 2> {
    public:
        using type = Tensor<Type>;
    };

    /*!\brief Not define, so do not use this class.*/
    template <class Type> class rankType<Type, 3> {
    public:
        using type = Tensor<Type>;
    };

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator+(const Tensor<Type> &t1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(t1.xx() + t2.xx(), t1.xy() + t2.xy(), t1.xz() + t2.xz(),
                            t1.yx() + t2.yx(), t1.yy() + t2.yy(), t1.yz() + t2.yz(),
                            t1.zx() + t2.zx(), t1.zy() + t2.zy(), t1.zz() + t2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator-(const Tensor<Type> &t2) {
        return Tensor<Type>(-t2.xx(), -t2.xy(), -t2.xz(), -t2.yx(), -t2.yy(), -t2.yz(), -t2.zx(),
                            -t2.zy(), -t2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator-(const Tensor<Type> &t1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(t1.xx() - t2.xx(), t1.xy() - t2.xy(), t1.xz() - t2.xz(),
                            t1.yx() - t2.yx(), t1.yy() - t2.yy(), t1.yz() - t2.yz(),
                            t1.zx() - t2.zx(), t1.zy() - t2.zy(), t1.zz() - t2.zz());
    }

    /**\brief Hodge Dual operator (tensor -> vector)*/
    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const Tensor<Type> &t) {
        return Vector<Type>(t.yz(), -t.xz(), t.xy());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator*(const Vector<Type> &v) {
        return Tensor<Type>(0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0);
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator*(const Tensor<Type> &t1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(t1.xx() * t2.xx() + t1.xy() * t2.yx() + t1.xz() * t2.zx(),
                            t1.xx() * t2.xy() + t1.xy() * t2.yy() + t1.xz() * t2.zy(),
                            t1.xx() * t2.xz() + t1.xy() * t2.yz() + t1.xz() * t2.zz(),

                            t1.yx() * t2.xx() + t1.yy() * t2.yx() + t1.yz() * t2.zx(),
                            t1.yx() * t2.xy() + t1.yy() * t2.yy() + t1.yz() * t2.zy(),
                            t1.yx() * t2.xz() + t1.yy() * t2.yz() + t1.yz() * t2.zz(),

                            t1.zx() * t2.xx() + t1.zy() * t2.yx() + t1.zz() * t2.zx(),
                            t1.zx() * t2.xy() + t1.zy() * t2.yy() + t1.zz() * t2.zy(),
                            t1.zx() * t2.xz() + t1.zy() * t2.yz() + t1.zz() * t2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator*(const Type t1, const Tensor<Type> &t2) {
        return Tensor<Type>(t1 * t2.xx(), t1 * t2.xy(), t1 * t2.xz(), t1 * t2.yx(), t1 * t2.yy(),
                            t1 * t2.yz(), t1 * t2.zx(), t1 * t2.zy(), t1 * t2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator&(const Tensor<Type> &t1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(t1.xx() * t2.xx() + t1.xy() * t2.yx() + t1.xz() * t2.zx(),
                            t1.xx() * t2.xy() + t1.xy() * t2.yy() + t1.xz() * t2.zy(),
                            t1.xx() * t2.xz() + t1.xy() * t2.yz() + t1.xz() * t2.zz(),

                            t1.yx() * t2.xx() + t1.yy() * t2.yx() + t1.yz() * t2.zx(),
                            t1.yx() * t2.xy() + t1.yy() * t2.yy() + t1.yz() * t2.zy(),
                            t1.yx() * t2.xz() + t1.yy() * t2.yz() + t1.yz() * t2.zz(),

                            t1.zx() * t2.xx() + t1.zy() * t2.yx() + t1.zz() * t2.zx(),
                            t1.zx() * t2.xy() + t1.zy() * t2.yy() + t1.zz() * t2.zy(),
                            t1.zx() * t2.xz() + t1.zy() * t2.yz() + t1.zz() * t2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const Tensor<Type> &t,
                                                          const Vector<Type> &v) {
        return Vector<Type>(t.xx() * v.x() + t.xy() * v.y() + t.xz() * v.z(),
                            t.yx() * v.x() + t.yy() * v.y() + t.yz() * v.z(),
                            t.zx() * v.x() + t.zy() * v.y() + t.zz() * v.z());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator&(const Vector<Type> &v1,
                                                          const Vector<Type> &v2) {
        return Tensor<Type>(v1.x() * v2.x(), v1.x() * v2.y(), v1.x() * v2.z(), v1.y() * v2.x(),
                            v1.y() * v2.y(), v1.y() * v2.z(), v1.z() * v2.x(), v1.z() * v2.y(),
                            v1.z() * v2.z());
    }

    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator&(const Tensor<Type> &t,
                                                          const Vector<Type> &v) {
        return Vector<Type>(t.xx() * v.x() + t.xy() * v.y() + t.xz() * v.z(),
                            t.yx() * v.x() + t.yy() * v.y() + t.yz() * v.z(),
                            t.zx() * v.x() + t.zy() * v.y() + t.zz() * v.z());
    }

    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> operator*(const Vector<Type> &v,
                                                          const Tensor<Type> &t) {
        return Vector<Type>(v.x() * t.xx() + v.y() * t.yx() + v.z() * t.zx(),
                            v.x() * t.xy() + v.y() * t.yy() + v.z() * t.zy(),
                            v.x() * t.xz() + v.y() * t.yz() + v.z() * t.zz());
    }

    template <class Type>
    inline Vector<Type> operator&(const Vector<Type> &v, const Tensor<Type> &t) {
        return Vector<Type>(v.x() * t.xx() + v.y() * t.yx() + v.z() * t.zx(),
                            v.x() * t.xy() + v.y() * t.yy() + v.z() * t.zy(),
                            v.x() * t.xz() + v.y() * t.yz() + v.z() * t.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator||(const Vector<Type> &v1,
                                                           const Vector<Type> &v2) {
        return Tensor<Type>(v1.x() * v2.x(), v1.x() * v2.y(), v1.x() * v2.z(), v1.y() * v2.x(),
                            v1.y() * v2.y(), v1.y() * v2.z(), v1.z() * v2.x(), v1.z() * v2.y(),
                            v1.z() * v2.z());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator||(const Tensor<Type> &v1,
                                                           const Tensor<Type> &v2) {
        return Tensor<Type>(v1.xx() * v2.xx(), v1.xy() * v2.xy(), v1.xz() * v2.xz(),
                            v1.yx() * v2.yx(), v1.yy() * v2.yy(), v1.yz() * v2.yz(),
                            v1.zx() * v2.zx(), v1.zy() * v2.zy(), v1.zz() * v2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline typename innerProduct<Vector<Type>, Tensor<Type>>::type
    operator/(const Vector<Type> &v, const Tensor<Type> &t) {
        return inv(t) * v;
    }

    /**\brief Return the trace of a tensor*/
    template <class Type> hur_nodiscard constexpr inline Type tr(const Tensor<Type> &t) {
        return t.xx() + t.yy() + t.zz();
    }

    /**\brief Return the spherical part of a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SphericalTensor<Type> sph(const Tensor<Type> &t) {
        return (1.0 / 3.0) * tr(t);
    }

    /**\brief Return the symmetric part of a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> symm(const Tensor<Type> &t) {
        return SymmTensor<Type>(t.xx(), 0.5 * (t.xy() + t.yx()), 0.5 * (t.xz() + t.zx()), t.yy(),
                                0.5 * (t.yz() + t.zy()), t.zz());
    }

    /**\brief Return twice the symmetric part of a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline SymmTensor<Type> twoSymm(const Tensor<Type> &t) {
        return SymmTensor<Type>(2 * t.xx(), (t.xy() + t.yx()), (t.xz() + t.zx()), 2 * t.yy(),
                                (t.yz() + t.zy()), 2 * t.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Vector<Type> diagToVector(const Tensor<Type> &t) {
        return Vector<Type>(t.xx(), t.yy(), t.zz());
    }

    /**\brief Return the skew-symmetric part of a tensor*/
    template <class Type> hur_nodiscard constexpr inline Tensor<Type> skew(const Tensor<Type> &t) {
        return Tensor<Type>(0.0, 0.5 * (t.xy() - t.yx()), 0.5 * (t.xz() - t.zx()),
                            0.5 * (t.yx() - t.xy()), 0.0, 0.5 * (t.yz() - t.zy()),
                            0.5 * (t.zx() - t.xz()), 0.5 * (t.zy() - t.yz()), 0.0);
    }

    /**\brief Return the magnitude of the skew-symmetric part of a tensor*/
    template <class Type> hur_nodiscard inline real skewMagnitude(const Tensor<Type> &t) {
        Type w1 = t.xy() - t.yx();
        Type w2 = t.xz() - t.zx();
        Type w3 = t.zy() - t.yz();
        return sqrt(real(w1 * w1 + w2 * w2 + w3 * w3));
    }

    /**\brief Return the square of the magnitude of the skew-symmetric part of a tensor*/
    template <class Type> hur_nodiscard constexpr inline real skewMagSqr(const Tensor<Type> &t) {
        Type w1 = t.xy() - t.yx();
        Type w2 = t.xz() - t.zx();
        Type w3 = t.zy() - t.yz();
        return real(w1 * w1 + w2 * w2 + w3 * w3);
    }

    /**\brief Return the skew-symmetric part of a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline const Tensor<Type> &skew(const SymmTensor<Type> &st) {
        return Tensor<Type>::zero;
    }

    /**\brief Return the magnitude of the skew-symmetric part of a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline real skewMagnitude(const SymmTensor<Type> &st) {
        return 0.0;
    }

    /**\brief Return the deviatoric part of a tensor*/
    template <class Type> hur_nodiscard inline Tensor<Type> dev(const Tensor<Type> &t) {
        return t - SphericalTensor<Type>::oneThirdI * tr(t);
    }

    /**\brief Return the deviatoric part of a tensor*/
    template <class Type> hur_nodiscard inline Tensor<Type> dev2(const Tensor<Type> &t) {
        return t - SphericalTensor<Type>::twoThirdsI * tr(t);
    }

    /**\brief Return the determinant of a tensor*/
    template <class Type> hur_nodiscard constexpr inline Type det(const Tensor<Type> &t) {
        return (t.xx() * t.yy() * t.zz() + t.xy() * t.yz() * t.zx() + t.xz() * t.yx() * t.zy() -
                t.xx() * t.yz() * t.zy() - t.xy() * t.yx() * t.zz() - t.xz() * t.yy() * t.zx());
    }

    /**\brief Return the cofactor tensor of a tensor*/
    template <class Type> hur_nodiscard constexpr inline Tensor<Type> cof(const Tensor<Type> &t) {
        return Tensor<Type>(t.yy() * t.zz() - t.zy() * t.yz(), t.zx() * t.yz() - t.yx() * t.zz(),
                            t.yx() * t.zy() - t.yy() * t.zx(),

                            t.xz() * t.zy() - t.xy() * t.zz(), t.xx() * t.zz() - t.xz() * t.zx(),
                            t.xy() * t.zx() - t.xx() * t.zy(),

                            t.xy() * t.yz() - t.xz() * t.yy(), t.yx() * t.xz() - t.xx() * t.yz(),
                            t.xx() * t.yy() - t.yx() * t.xy());
    }

    /**\brief Return the inverse of a tensor given the determinant*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> inv(const Tensor<Type> &t, const Type dett) {
        return Tensor<Type>(t.yy() * t.zz() - t.zy() * t.yz(), t.xz() * t.zy() - t.xy() * t.zz(),
                            t.xy() * t.yz() - t.xz() * t.yy(),

                            t.zx() * t.yz() - t.yx() * t.zz(), t.xx() * t.zz() - t.xz() * t.zx(),
                            t.yx() * t.xz() - t.xx() * t.yz(),

                            t.yx() * t.zy() - t.yy() * t.zx(), t.xy() * t.zx() - t.xx() * t.zy(),
                            t.xx() * t.yy() - t.yx() * t.xy()) /
               dett;
    }
    template <class Type> hur_nodiscard constexpr inline Tensor<Type> inv(const Tensor<Type> &t) {
        return inv(t, det(t));
    }

    template <class Type> hur_nodiscard constexpr inline Type invariantI(const Tensor<Type> &t) {
        return tr(t);
    }

    /**\brief Return the 2nd invariant of a tensor*/
    template <class Type> hur_nodiscard constexpr inline Type invariantII(const Tensor<Type> &t) {
        return (t.xx() * t.yy() + t.yy() * t.zz() + t.xx() * t.zz() - t.xy() * t.yx() -
                t.yz() * t.zy() - t.xz() * t.zx());
    }

    /**\brief Return the 3rd invariant of a tensor*/
    template <class Type> hur_nodiscard constexpr inline Type invariantIII(const Tensor<Type> &t) {
        return det(t);
    }

    /**\brief sum aibi between two tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aibi(const Tensor<Type> &t1, const Tensor<Type> &t2) {
        return Type(t1.xx() * t2.xx() + t1.yy() * t2.yy() + t1.zz() * t2.zz());
    }

    /**\brief sum aibi between a symmetric tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aibi(const SymmTensor<Type> &st1, const Tensor<Type> &t2) {
        return Type(st1.xx() * t2.xx() + st1.yy() * t2.yy() + st1.zz() * t2.zz());
    }

    /**\brief sum aibi between a tensor and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aibi(const Tensor<Type> &t1, const SymmTensor<Type> &st2) {
        return Type(t1.xx() * st2.xx() + t1.yy() * st2.yy() + t1.zz() * st2.zz());
    }

    /**\brief sum aibi between two symmetric tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aibi(const SymmTensor<Type> &t1,
                                             const SymmTensor<Type> &st2) {
        return Type(t1.xx() * st2.xx() + t1.yy() * st2.yy() + t1.zz() * st2.zz());
    }

    /**\brief sum aibi between a spherical tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aibi(const SphericalTensor<Type> &st1,
                                             const Tensor<Type> &t2) {
        return (st1.ii() * t2.xx() + st1.ii() * t2.yy() + st1.ii() * t2.zz());
    }

    /**\brief sum aibi between a tensor and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aibi(const Tensor<Type> &t1,
                                             const SphericalTensor<Type> &st2) {
        return Type(t1.xx() * st2.ii() + t1.yy() * st2.ii() + t1.zz() * st2.ii());
    }

    /**\brief sum aibi between two spherical tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aibi(const SphericalTensor<Type> &t1,
                                             const SphericalTensor<Type> &st2) {
        return Type(t1.ii() * st2.ii() + t1.ii() * st2.ii() + t1.ii() * st2.ii());
    }

    /**\brief sum aibi between a diagnal tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aibi(const DiagTensor<Type> &st1, const Tensor<Type> &t2) {
        return Type(st1.xx() * t2.xx() + st1.yy() * t2.yy() + st1.zz() * t2.zz());
    }

    /**\brief sum aibi between a tensor and a diagnal tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aibi(const Tensor<Type> &t1, const DiagTensor<Type> &st2) {
        return Type(t1.xx() * st2.xx() + t1.yy() * st2.yy() + t1.zz() * st2.zz());
    }

    /**\brief sum aibi between two diagnal tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aibi(const DiagTensor<Type> &t1,
                                             const DiagTensor<Type> &st2) {
        return Type(t1.xx() * st2.xx() + t1.yy() * st2.yy() + t1.zz() * st2.zz());
    }

    // aijbi-----------------------------

    /**\brief sum aijbi between two tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbi(const Tensor<Type> &t1, const Tensor<Type> &t2) {
        return Type((t1.xx() + t1.yx() + t1.zx()) * t2.xx() +
                    (t1.xy() + t1.yy() + t1.zy()) * t2.yy() +
                    (t1.xz() + t1.yz() + t1.zz()) * t2.zz());
    }

    /**\brief sum aijbi between a symmetric tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbi(const SymmTensor<Type> &st1, const Tensor<Type> &t2) {
        return Type((st1.xx() + st1.xy() + st1.xz()) * t2.xx() +
                    (st1.xy() + st1.yy() + st1.yz()) * t2.yy() +
                    (st1.xz() + st1.yz() + st1.zz()) * t2.zz());
    }

    /**\brief sum aijbi between a tensor and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbi(const Tensor<Type> &t1, const SymmTensor<Type> &t2) {
        return Type((t1.xx() + t1.yx() + t1.zx()) * t2.xx() +
                    (t1.xy() + t1.yy() + t1.zy()) * t2.yy() +
                    (t1.xz() + t1.yz() + t1.zz()) * t2.zz());
    }

    /**\brief sum aijbi between two symmetric tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbi(const SymmTensor<Type> &st1,
                                              const SymmTensor<Type> &st2) {
        return Type((st1.xx() + st1.xy() + st1.xz()) * st2.xx() +
                    (st1.xy() + st1.yy() + st1.yz()) * st2.yy() +
                    (st1.xz() + st1.yz() + st1.zz()) * st2.zz());
    }

    /**\brief sum aijbi between a spherical tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbi(const SphericalTensor<Type> &st1,
                                              const Tensor<Type> &t2) {
        return Type(st1.ii() * t2.xx() + st1.ii() * t2.yy() + st1.ii() * t2.zz());
    }

    /**\brief sum aijbi between a tensor and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbi(const Tensor<Type> &t1,
                                              const SphericalTensor<Type> &st2) {
        return Type(t1.xx() * st2.ii() + t1.yy() * st2.ii() + t1.zz() * st2.ii());
    }

    /**\brief sum aijbi between two spherical tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbi(const SphericalTensor<Type> &t1,
                                              const SphericalTensor<Type> &st2) {
        return Type(t1.ii() * st2.ii() + t1.ii() * st2.ii() + t1.ii() * st2.ii());
    }

    /**\brief sum aijbi between a diagnal tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbi(const DiagTensor<Type> &st1, const Tensor<Type> &t2) {
        return Type(st1.xx() * t2.xx() + st1.yy() * t2.yy() + st1.zz() * t2.zz());
    }

    /**\brief sum aijbi between a tensor and a diagnal tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbi(const Tensor<Type> &t1, const DiagTensor<Type> &st2) {
        return Type(t1.xx() * st2.xx() + t1.yy() * st2.yy() + t1.zz() * st2.zz());
    }

    /**\brief sum aijbi between two diagnal tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbi(const DiagTensor<Type> &t1,
                                              const DiagTensor<Type> &st2) {
        return Type(t1.xx() * st2.xx() + t1.yy() * st2.yy() + t1.zz() * st2.zz());
    }

    // aii2-------------------------------------------

    /**\brief sum aii2 for tensor*/
    template <class Type> hur_nodiscard constexpr inline Type aii2(const Tensor<Type> &t1) {
        return Type(sqr(t1.xx()) + sqr(t1.yy()) + sqr(t1.zz()));
    }

    /**\brief sum aii2 for spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aii2(const SphericalTensor<Type> &t1) {
        return Type(sqr(t1.ii()) + sqr(t1.ii()) + sqr(t1.ii()));
    }

    /**\brief sum aii2 for symmetric tensor*/
    template <class Type> hur_nodiscard constexpr inline Type aii2(const SymmTensor<Type> &t1) {
        return Type(sqr(t1.xx()) + sqr(t1.yy()) + sqr(t1.zz()));
    }

    /**\brief sum aii2 for diagnal tensor*/
    template <class Type> hur_nodiscard constexpr inline Type aii2(const DiagTensor<Type> &t1) {
        return Type(sqr(t1.xx()) + sqr(t1.yy()) + sqr(t1.zz()));
    }

    // aijbij--------------------------------

    /**\brief sum aijbij between two tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbij(const Tensor<Type> &t1, const Tensor<Type> &t2) {
        return Type(t1.xx() * t2.xx() + t1.xy() * t2.xy() + t1.xz() * t2.xz() + t1.yx() * t2.yx() +
                    t1.yy() * t2.yy() + t1.yz() * t2.yz() + t1.zx() * t2.zx() + t1.zy() * t2.zy() +
                    t1.zz() * t2.zz());
    }

    /**\brief sum aijbij between a symmetric tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbij(const SymmTensor<Type> &st1,
                                               const Tensor<Type> &t2) {
        return Type(st1.xx() * t2.xx() + st1.xy() * t2.xy() + st1.xz() * t2.xz() +
                    st1.xy() * t2.yx() + st1.yy() * t2.yy() + st1.yz() * t2.yz() +
                    st1.xz() * t2.zx() + st1.yz() * t2.zy() + st1.zz() * t2.zz());
    }

    /**\brief sum aijbij between a tensor and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbij(const Tensor<Type> &t1,
                                               const SymmTensor<Type> &st2) {
        return Type(t1.xx() * st2.xx() + t1.xy() * st2.xy() + t1.xz() * st2.xz() +
                    t1.yx() * st2.xy() + t1.yy() * st2.yy() + t1.yz() * st2.yz() +
                    t1.zx() * st2.xz() + t1.zy() * st2.yz() + t1.zz() * st2.zz());
    }

    /**\brief sum aijbij between two symmetric tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbij(const SymmTensor<Type> &st1,
                                               const SymmTensor<Type> &st2) {
        return Type(st1.xx() * st2.xx() + 2.0 * st1.xy() * st2.xy() + 2.0 * st1.xz() * st2.xz() +
                    st1.yy() * st2.yy() + 2.0 * st1.yz() * st2.yz() + st1.zz() * st2.zz());
    }

    /**\brief sum aijbij between a spherical tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbij(const SphericalTensor<Type> &st1,
                                               const Tensor<Type> &t2) {
        return Type(st1.ii() * t2.xx() + st1.ii() * t2.yy() + st1.ii() * t2.zz());
    }

    /**\brief sum aijbij between a tensor and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbij(const Tensor<Type> &t1,
                                               const SphericalTensor<Type> &st2) {
        return Type(t1.xx() * st2.ii() + t1.yy() * st2.ii() + t1.zz() * st2.ii());
    }

    /**\brief sum aijbij between two spherical tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbij(const SphericalTensor<Type> &t1,
                                               const SphericalTensor<Type> &st2) {
        return Type(t1.ii() * st2.ii() + t1.ii() * st2.ii() + t1.ii() * st2.ii());
    }

    /**\brief sum aijbij between a diagnal tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbij(const DiagTensor<Type> &st1,
                                               const Tensor<Type> &t2) {
        return Type(st1.xx() * t2.xx() + st1.yy() * t2.yy() + st1.zz() * t2.zz());
    }

    /**\brief sum aijbij between a tensor and a diagnal tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbij(const Tensor<Type> &t1,
                                               const DiagTensor<Type> &st2) {
        return Type(t1.xx() * st2.xx() + t1.yy() * st2.yy() + t1.zz() * st2.zz());
    }

    /**\brief sum aijbij between two diagnal tensors*/
    template <class Type>
    hur_nodiscard constexpr inline Type aijbij(const DiagTensor<Type> &t1,
                                               const DiagTensor<Type> &st2) {
        return Type(t1.xx() * st2.xx() + t1.yy() * st2.yy() + t1.zz() * st2.zz());
    }

    template <class Type>
    inline Tensor<Type> hur_nodiscard constexpr operator+(const SphericalTensor<Type> &st1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(st1.ii() + t2.xx(), t2.xy(), t2.xz(), t2.yx(), st1.ii() + t2.yy(),
                            t2.yz(), t2.zx(), t2.zy(), st1.ii() + t2.zz());
    }

    template <class Type>
    inline Tensor<Type> hur_nodiscard constexpr operator+(const Tensor<Type> &t1,
                                                          const SphericalTensor<Type> &st2) {
        return Tensor<Type>(t1.xx() + st2.ii(), t1.xy(), t1.xz(), t1.yx(), t1.yy() + st2.ii(),
                            t1.yz(), t1.zx(), t1.zy(), t1.zz() + st2.ii());
    }

    template <class Type>
    inline Tensor<Type> hur_nodiscard constexpr operator-(const SphericalTensor<Type> &st1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(st1.ii() - t2.xx(), -t2.xy(), -t2.xz(), -t2.yx(), st1.ii() - t2.yy(),
                            -t2.yz(), -t2.zx(), -t2.zy(), st1.ii() - t2.zz());
    }

    template <class Type>
    inline Tensor<Type> hur_nodiscard constexpr operator-(const Tensor<Type> &t1,
                                                          const SphericalTensor<Type> &st2) {
        return Tensor<Type>(t1.xx() - st2.ii(), t1.xy(), t1.xz(), t1.yx(), t1.yy() - st2.ii(),
                            t1.yz(), t1.zx(), t1.zy(), t1.zz() - st2.ii());
    }

    /**\brief Inner-product between a spherical tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator&(const SphericalTensor<Type> &st1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(st1.ii() * t2.xx(), st1.ii() * t2.xy(), st1.ii() * t2.xz(),
                            st1.ii() * t2.yx(), st1.ii() * t2.yy(), st1.ii() * t2.yz(),
                            st1.ii() * t2.zx(), st1.ii() * t2.zy(), st1.ii() * t2.zz());
    }

    /**\brief Inner-product between a spherical tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator*(const SphericalTensor<Type> &st1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(st1.ii() * t2.xx(), st1.ii() * t2.xy(), st1.ii() * t2.xz(),
                            st1.ii() * t2.yx(), st1.ii() * t2.yy(), st1.ii() * t2.yz(),
                            st1.ii() * t2.zx(), st1.ii() * t2.zy(), st1.ii() * t2.zz());
    }

    /**\brief Inner-product between a tensor and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator&(const Tensor<Type> &t1,
                                                          const SphericalTensor<Type> &st2) {
        return Tensor<Type>(t1.xx() * st2.ii(), t1.xy() * st2.ii(), t1.xz() * st2.ii(),
                            t1.yx() * st2.ii(), t1.yy() * st2.ii(), t1.yz() * st2.ii(),
                            t1.zx() * st2.ii(), t1.zy() * st2.ii(), t1.zz() * st2.ii());
    }

    /**\brief Inner-product between a tensor and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator*(const Tensor<Type> &t1,
                                                          const SphericalTensor<Type> &st2) {
        return Tensor<Type>(t1.xx() * st2.ii(), t1.xy() * st2.ii(), t1.xz() * st2.ii(),
                            t1.yx() * st2.ii(), t1.yy() * st2.ii(), t1.yz() * st2.ii(),
                            t1.zx() * st2.ii(), t1.zy() * st2.ii(), t1.zz() * st2.ii());
    }

    /**\brief Double-dot-product between a spherical tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type operator&&(const SphericalTensor<Type> &st1,
                                                   const Tensor<Type> &t2) {
        return (st1.ii() * t2.xx() + st1.ii() * t2.yy() + st1.ii() * t2.zz());
    }

    /**\brief Double-dot-product between a tensor and a spherical tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type operator&&(const Tensor<Type> &t1,
                                                   const SphericalTensor<Type> &st2) {
        return (t1.xx() * st2.ii() + t1.yy() * st2.ii() + t1.zz() * st2.ii());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator+(const SymmTensor<Type> &st1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(st1.xx() + t2.xx(), st1.xy() + t2.xy(), st1.xz() + t2.xz(),
                            st1.xy() + t2.yx(), st1.yy() + t2.yy(), st1.yz() + t2.yz(),
                            st1.xz() + t2.zx(), st1.yz() + t2.zy(), st1.zz() + t2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator+(const Tensor<Type> &t1,
                                                          const SymmTensor<Type> &st2) {
        return Tensor<Type>(t1.xx() + st2.xx(), t1.xy() + st2.xy(), t1.xz() + st2.xz(),
                            t1.yx() + st2.xy(), t1.yy() + st2.yy(), t1.yz() + st2.yz(),
                            t1.zx() + st2.xz(), t1.zy() + st2.yz(), t1.zz() + st2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator-(const SymmTensor<Type> &st1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(st1.xx() - t2.xx(), st1.xy() - t2.xy(), st1.xz() - t2.xz(),
                            st1.xy() - t2.yx(), st1.yy() - t2.yy(), st1.yz() - t2.yz(),
                            st1.xz() - t2.zx(), st1.yz() - t2.zy(), st1.zz() - t2.zz());
    }

    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator-(const Tensor<Type> &t1,
                                                          const SymmTensor<Type> &st2) {
        return Tensor<Type>(t1.xx() - st2.xx(), t1.xy() - st2.xy(), t1.xz() - st2.xz(),
                            t1.yx() - st2.xy(), t1.yy() - st2.yy(), t1.yz() - st2.yz(),
                            t1.zx() - st2.xz(), t1.zy() - st2.yz(), t1.zz() - st2.zz());
    }

    /**\brief Inner-product between a symmetric tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator&(const SymmTensor<Type> &st1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(st1.xx() * t2.xx() + st1.xy() * t2.yx() + st1.xz() * t2.zx(),
                            st1.xx() * t2.xy() + st1.xy() * t2.yy() + st1.xz() * t2.zy(),
                            st1.xx() * t2.xz() + st1.xy() * t2.yz() + st1.xz() * t2.zz(),

                            st1.xy() * t2.xx() + st1.yy() * t2.yx() + st1.yz() * t2.zx(),
                            st1.xy() * t2.xy() + st1.yy() * t2.yy() + st1.yz() * t2.zy(),
                            st1.xy() * t2.xz() + st1.yy() * t2.yz() + st1.yz() * t2.zz(),

                            st1.xz() * t2.xx() + st1.yz() * t2.yx() + st1.zz() * t2.zx(),
                            st1.xz() * t2.xy() + st1.yz() * t2.yy() + st1.zz() * t2.zy(),
                            st1.xz() * t2.xz() + st1.yz() * t2.yz() + st1.zz() * t2.zz());
    }
    /**\brief Inner-product between a symmetric tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator*(const SymmTensor<Type> &st1,
                                                          const Tensor<Type> &t2) {
        return Tensor<Type>(st1.xx() * t2.xx() + st1.xy() * t2.yx() + st1.xz() * t2.zx(),
                            st1.xx() * t2.xy() + st1.xy() * t2.yy() + st1.xz() * t2.zy(),
                            st1.xx() * t2.xz() + st1.xy() * t2.yz() + st1.xz() * t2.zz(),

                            st1.xy() * t2.xx() + st1.yy() * t2.yx() + st1.yz() * t2.zx(),
                            st1.xy() * t2.xy() + st1.yy() * t2.yy() + st1.yz() * t2.zy(),
                            st1.xy() * t2.xz() + st1.yy() * t2.yz() + st1.yz() * t2.zz(),

                            st1.xz() * t2.xx() + st1.yz() * t2.yx() + st1.zz() * t2.zx(),
                            st1.xz() * t2.xy() + st1.yz() * t2.yy() + st1.zz() * t2.zy(),
                            st1.xz() * t2.xz() + st1.yz() * t2.yz() + st1.zz() * t2.zz());
    }

    /**\brief Inner-product between a tensor and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator&(const Tensor<Type> &t1,
                                                          const SymmTensor<Type> &st2) {
        return Tensor<Type>(t1.xx() * st2.xx() + t1.xy() * st2.xy() + t1.xz() * st2.xz(),
                            t1.xx() * st2.xy() + t1.xy() * st2.yy() + t1.xz() * st2.yz(),
                            t1.xx() * st2.xz() + t1.xy() * st2.yz() + t1.xz() * st2.zz(),

                            t1.yx() * st2.xx() + t1.yy() * st2.xy() + t1.yz() * st2.xz(),
                            t1.yx() * st2.xy() + t1.yy() * st2.yy() + t1.yz() * st2.yz(),
                            t1.yx() * st2.xz() + t1.yy() * st2.yz() + t1.yz() * st2.zz(),

                            t1.zx() * st2.xx() + t1.zy() * st2.xy() + t1.zz() * st2.xz(),
                            t1.zx() * st2.xy() + t1.zy() * st2.yy() + t1.zz() * st2.yz(),
                            t1.zx() * st2.xz() + t1.zy() * st2.yz() + t1.zz() * st2.zz());
    }
    /**\brief Inner-product between a tensor and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator*(const Tensor<Type> &t1,
                                                          const SymmTensor<Type> &st2) {
        return Tensor<Type>(t1.xx() * st2.xx() + t1.xy() * st2.xy() + t1.xz() * st2.xz(),
                            t1.xx() * st2.xy() + t1.xy() * st2.yy() + t1.xz() * st2.yz(),
                            t1.xx() * st2.xz() + t1.xy() * st2.yz() + t1.xz() * st2.zz(),

                            t1.yx() * st2.xx() + t1.yy() * st2.xy() + t1.yz() * st2.xz(),
                            t1.yx() * st2.xy() + t1.yy() * st2.yy() + t1.yz() * st2.yz(),
                            t1.yx() * st2.xz() + t1.yy() * st2.yz() + t1.yz() * st2.zz(),

                            t1.zx() * st2.xx() + t1.zy() * st2.xy() + t1.zz() * st2.xz(),
                            t1.zx() * st2.xy() + t1.zy() * st2.yy() + t1.zz() * st2.yz(),
                            t1.zx() * st2.xz() + t1.zy() * st2.yz() + t1.zz() * st2.zz());
    }

    /**\brief Double-dot-product between a symmetric tensor and a tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type operator&&(const SymmTensor<Type> &st1,
                                                   const Tensor<Type> &t2) {
        return (st1.xx() * t2.xx() + st1.xy() * t2.xy() + st1.xz() * t2.xz() + st1.xy() * t2.yx() +
                st1.yy() * t2.yy() + st1.yz() * t2.yz() + st1.xz() * t2.zx() + st1.yz() * t2.zy() +
                st1.zz() * t2.zz());
    }

    /**\brief Double-dot-product between a tensor and a symmetric tensor*/
    template <class Type>
    hur_nodiscard constexpr inline Type operator&&(const Tensor<Type> &t1,
                                                   const SymmTensor<Type> &st2) {
        return (t1.xx() * st2.xx() + t1.xy() * st2.xy() + t1.xz() * st2.xz() + t1.yx() * st2.xy() +
                t1.yy() * st2.yy() + t1.yz() * st2.yz() + t1.zx() * st2.xz() + t1.zy() * st2.yz() +
                t1.zz() * st2.zz());
    }

    /**\brief kI + A
            A is a tensor; I is a identity tensor; k is a real*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator+(const Type &k, const Tensor<Type> &t2) {
        return Tensor<Type>(k + t2.xx(), t2.xy(), t2.xz(), t2.yx(), k + t2.yy(), t2.yz(), t2.zx(),
                            t2.zy(), k + t2.zz());
    }

    /**\brief A + kI
            A is a tensor; I is a identity tensor; k is a real*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator+(const Tensor<Type> &t1, const Type &k) {
        return Tensor<Type>(t1.xx() + k, t1.xy(), t1.xz(), t1.yx(), t1.yy() + k, t1.yz(), t1.zx(),
                            t1.zy(), t1.zz() + k);
    }

    /**\brief kI - A
            A is a tensor; I is a identity tensor; k is a real*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator-(const Type &k, const Tensor<Type> &t2) {
        return Tensor<Type>(k - t2.xx(), -t2.xy(), -t2.xz(), -t2.yx(), k - t2.yy(), -t2.yz(),
                            -t2.zx(), -t2.zy(), k - t2.zz());
    }

    /**\brief A - kI
            A is a tensor; I is a identity tensor; k is a real*/
    template <class Type>
    hur_nodiscard constexpr inline Tensor<Type> operator-(const Tensor<Type> &t1, const Type &k) {
        return Tensor<Type>(t1.xx() - k, t1.xy(), t1.xz(), t1.yx(), t1.yy() - k, t1.yz(), t1.zx(),
                            t1.zy(), t1.zz() - k);
    }
} //  namespace OpenHurricane
