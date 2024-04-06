/*!
 * \file Array.hpp
 * \brief Headers of the Array.
 *        The subroutines and functions are in the <i>Array.inl</i> file.
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

#include "List.hpp"
#include "dataStructure.hpp"
#include "fileIsstream.hpp"
#include "fileOsstream.hpp"

namespace OpenHurricane {

    /*!\brief The template class of Array.*/
    template <class Type> class Array : public List<Type> {
    public:
        using Base = List<Type>;
        using value_type = typename Base::value_type;
        using reference = typename Base::reference;
        using const_reference = typename Base::const_reference;
        using difference_type = typename Base::difference_type;
        using size_type = typename List<Type>::size_type;

        /*!\brief Component type.*/
        using elementType = typename feature<Type>::elementType;

        /*!\brief Declare type of subArray.*/
        using subArray = Array<Type>;

    protected:
        inline virtual void setCurRhs(const Type &curRhs) const {}

    public:
        static constexpr size_type npos = Base::npos;

        hur_nodiscard inline static const Array<Type> &nullObject() {
            return NullRefObj::nullRef<Array<Type>>();
        }

        inline Array() : Base() {}
        inline explicit Array(const size_type s) : Base(s) {}

        inline Array(const size_type s, const Type &ele) : Base(s, ele) {}
        inline Array &operator=(const Type &ele) {
            Base::operator=(ele);
            return *this;
        }

        inline Array(const size_type s, const zero) : Base(s, Zero) {}
        inline Array &operator=(const zero) {
            Base::operator=(Zero);
            return *this;
        }

        inline Array(const Array &other) : Base(other) {}
        inline Array &operator=(const Array &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        inline Array(Array &&other) noexcept : Base(std::move(other)) {}
        inline Array &operator=(Array &&other) noexcept {
            Base::operator=(std::move(other));
            return *this;
        }

        inline explicit Array(const Base &other) : Base(other) {}
        inline Array &operator=(const Base &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        hur_nodiscard inline uniquePtr<Array> clone() const {
            return uniquePtr<Array>(new Array(*this));
        }

        /*!\brief Return a sub-array.*/
        hur_nodiscard inline Array sub(const size_type Off = 0, const size_type Count = npos) {
            return Array(Base::sub(Off, Count));
        }

        /*!\brief Destructor.*/
        inline virtual ~Array() noexcept {}

        /*!
         * \brief Replace a component Array of the Array
         *  Only for Type which has components such as vector.
         */
        void replace(const int i, const Array<elementType> &l) {
            for (size_type j = 0; j < this->size(); ++j) {
                this->operator[](j)[i] = l[j];
            }
        }

        /*!
         * \brief Replace a component Array of the Array
         *  Only for Type which has components such as vector.
         */
        void replace(const int i, const elementType &c) {
            for (size_type j = 0; j < this->size(); j++) {
                this->operator[](j)[i] = c;
            }
        }

        /*!
         * \brief Set a component Array of the Array with the given component value.
         *  Only for Type which has components such as vector.
         */
        void setComponent(const int d, const elementType &c) {
            for (size_type j = 0; j < this->size(); j++) {
                this->operator[](j)[d] = c;
            }
        }

        /*!
         * \brief Set a component Array of the Array with the given component value.
         *  Only for Type which has components such as vector.
         */
        void setComponent(const int d, const Array<elementType> &c) {
            for (size_type j = 0; j < min(this->size(), c.size()); j++) {
                this->operator[](j)[d] = c[j];
            }
        }

        /*!
         * \brief Set all component Array of the Array with the given component value.
         *  Only for Type which has components such as vector.
         */
        void setComponent(const elementType &c) {
            for (size_type j = 0; j < this->size(); j++) {
                for (int i = 0; i < feature<Type>::nElements_; ++i) {
                    this->operator[](j)[i] = c;
                }
            }
        }

        hur_nodiscard Type average() const {
            if (this->size() == 0) {
                LFatal("Attempt to average null Array");
            }
            Type s = Zero;
            for (size_type i = 0; i < this->size(); ++i) {
                s += this->operator[](i);
            }
            s /= real(this->size());
            return s;
        }

        hur_nodiscard Type weightedAverage(const Array<real> &wF) const {
            if (this->size() == 0) {
                LFatal("Attempt to average null Array");
            }
            Type s = Zero;
            real w = Zero;
            for (size_type i = 0; i < this->size(); ++i) {
                s += wF[i] * this->operator[](i);
                w += wF[i];
            }
            s /= w;
            return s;
        }

        inline void clear() noexcept { Base::clear(); }

        hur_nodiscard Array<typename Array<Type>::elementType> component(const int d) const {
            Array<typename Array<Type>::elementType> comp(this->size());
            for (size_type i = 0; i < this->size(); ++i) {
                comp[i] = this->operator[](i)[d];
            }
            return comp;
        }

        inline void writeAveToPout(fileOsstream &fos, const Array<Type> &rhs, const Array<real> &cV,
                                   const size_type n, const size_type allN, const Type &rhs0,
                                   const bool calRhs0, const bool modifyRhs0 = false) const {}

        OpenHurricane::fileOsstream &writeToStream(fileOsstream &fos) const {
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                for (size_type i = 0; i < this->size(); ++i) {
                    fos.os() << this->operator[](i) << " ";
                }
            } else {
                if (this->size()) {
                    fos.write(reinterpret_cast<const char *>(&this->operator[](0)),
                              this->byteSize());
                }
            }
            return fos;
        }

        /*!\brief Start from v_[0].*/
        OpenHurricane::fileOsstream &writeToStream(fileOsstream &fos, const size_type outSize) const {
            size_type minSize = min(outSize, this->size());
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                for (size_type i = 0; i < minSize; ++i) {
                    fos.os() << this->operator[](i) << " ";
                }
            } else {
                if (minSize) {
                    fos.write(reinterpret_cast<const char *>(&this->operator[](0)),
                              minSize * sizeof(Type));
                }
            }
            return fos;
        }

        OpenHurricane::fileOsstream &writeToStreamWithFactor(fileOsstream &fos) const {
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                for (size_type i = 0; i < this->size(); ++i) {
                    fos.os() << this->operator[](i) << " ";
                }
            } else {
                if (this->size()) {
                    Array<Type> newF;
                    newF = this->clone();
                    fos.write(reinterpret_cast<const char *>(&newF[0]), this->byteSize());
                }
            }
            return fos;
        }

        /*!\brief Start from v_[0].*/
        OpenHurricane::fileOsstream &writeToStreamWithFactor(fileOsstream &fos,
                                                         const size_type outSize) const {
            size_type minSize = min(outSize, this->size());
            if (fos.format() == IOsstream::ASCII_FORMAT) {
                for (size_type i = 0; i < minSize; ++i) {
                    fos.os() << this->operator[](i) << " ";
                }
            } else {
                if (minSize) {
                    Array<Type> newF;
                    newF == this->clone();
                    fos.write(reinterpret_cast<const char *>(&newF[0]), minSize * sizeof(Type));
                }
            }
            return fos;
        }

        /*!\brief Start from v_[0].*/
        inline void writeMinMaxToStream(fileOsstream &fos, const size_type findSize) const {}
        inline void writeMinMaxToStreamWithFactor(fileOsstream &fos,
                                                  const size_type findSize) const {}
        inline void writeMinMaxToStreamWithFactorByMaster(fileOsstream &fos,
                                                          const size_type findSize) const {}

        void operator+=(const Array &f) {
            if (this->size() != f.size()) {
                LFatal("Attempt to add two Array in different size: %d != %d", this->size(),
                       f.size());
            }
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) += f[i];
            }
        }

        void operator+=(const List<real> &f) {
            if (this->size() != f.size()) {
                LFatal("Attempt to add two Array in different size: %d != %d", this->size(),
                       f.size());
            }
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) += f[i];
            }
        }

        void operator+=(const Type &t) {
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) += t;
            }
        }

        void operator-=(const Array &f) {
            if (this->size() != f.size()) {
                LFatal("Attempt to subtract two Array in different size: %d != %d", this->size(),
                       f.size());
            }
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) -= f[i];
            }
        }

        void operator-=(const List<real> &f) {
            if (this->size() != f.size()) {
                LFatal("Attempt to subtract two Array in different size: %d != %d", this->size(),
                       f.size());
            }
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) -= f[i];
            }
        }

        void operator-=(const Type &t) {
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) -= t;
            }
        }

        void operator*=(const Array &f) {
            if (this->size() != f.size()) {
                LFatal("Attempt to multiply two Array in different size: %d != %d", this->size(),
                       f.size());
            }
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) *= f[i];
            }
        }

        void operator*=(const List<real> &f) {
            if (this->size() != f.size()) {
                LFatal("Attempt to multiply two Array in different size: %d != %d", this->size(),
                       f.size());
            }
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) *= f[i];
            }
        }

        void operator*=(const real &t) {
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) *= t;
            }
        }

        void operator/=(const Array &f) {
            if (this->size() != f.size()) {
                LFatal("Attempt to devide two Array in different size: %d != %d", this->size(),
                       f.size());
            }
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) /= f[i];
            }
        }

        void operator/=(const List<real> &f) {
            if (this->size() != f.size()) {
                LFatal("Attempt to devide two Array in different size: %d != %d", this->size(),
                       f.size());
            }
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) /= f[i];
            }
        }

        void operator/=(const real &t) {
            if (t == real(0)) {
                LFatal("Attempt to devide Array by zero");
            }
            for (size_type i = 0; i < this->size(); ++i) {
                this->operator[](i) /= t;
            }
        }
    };

    template <class Type>
    hur_nodiscard inline Array<Type> operator+(const Array<Type> &f1, const Array<Type> &f2) {
        checkArraysSize(f1, f2, "+");
        Array<Type> f(f1.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] + f2[i];
        }
        return f;
    }
    template <class Type>
    hur_nodiscard inline Array<Type> operator+(Array<Type> &&f1, const Array<Type> &f2) noexcept {
        checkArraysSize(f1, f2, "+");
        Array<Type> tf(std::move(f1));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] += f2[i];
        }
        return tf;
    }
    template <class Type>
    hur_nodiscard inline Array<Type> operator+(const Array<Type> &f1, Array<Type> &&f2) noexcept {
        checkArraysSize(f1, f2, "+");
        Array<Type> tf(std::move(f2));
        for (typename Array<Type>::size_type i = 0; i < f1.size(); ++i) {
            tf[i] = f1[i] + tf[i];
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator+(Array<Type> &&f1, Array<Type> &&f2) noexcept {
        checkArraysSize(f1, f2, "+");
        Array<Type> tf(std::move(f1));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] += f2[i];
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator+(const Array<Type> &f, const Type &t) {
        Array<Type> ft(f.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            ft[i] = f[i] + t;
        }

        return ft;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator+(const Type &t, const Array<Type> &f) {
        return (f + t);
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator+(Array<Type> &&f, const Type &t) noexcept {
        Array<Type> tf(std::move(f));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] += t;
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator+(const Type &t, Array<Type> &&f) noexcept {
        Array<Type> tf(std::move(f));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] += t;
        }
        return tf;
    }

    template <class Type> hur_nodiscard inline Array<Type> operator-(const Array<Type> &f) {
        Array<Type> negF(f.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            negF[i] = -f[i];
        }
        return negF;
    }

    template <class Type> hur_nodiscard inline Array<Type> operator-(Array<Type> &&f) noexcept {
        Array<Type> tf(std::move(f));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = -tf[i];
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator-(const Array<Type> &f1, const Array<Type> &f2) {
        checkArraysSize(f1, f2, "-");
        Array<Type> f(f1.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] - f2[i];
        }

        return f;
    }
    template <class Type>
    hur_nodiscard inline Array<Type> operator-(Array<Type> &&f1, const Array<Type> &f2) noexcept {
        checkArraysSize(f1, f2, "-");
        Array<Type> tf(std::move(f1));
        for (integer i = 0; i < tf.size(); ++i) {
            tf[i] -= f2[i];
        }

        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator-(const Array<Type> &f1, Array<Type> &&f2) noexcept {
        checkArraysSize(f1, f2, "-");
        Array<Type> tf(std::move(f2));
        for (typename Array<Type>::size_type i = 0; i < f1.size(); ++i) {
            tf[i] = f1[i] - tf[i];
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator-(Array<Type> &&f1, Array<Type> &&f2) noexcept {
        checkArraysSize(f1, f2, "-");
        Array<Type> tf(std::move(f1));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] -= f2[i];
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator-(const Array<Type> &f, const Type &t) {
        Array<Type> ft(f.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            ft[i] = f[i] - t;
        }
        return ft;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator-(const Type &t, const Array<Type> &f) {
        Array<Type> tf(f.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            tf[i] = t - f[i];
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator-(Array<Type> &&f, const Type &t) noexcept {
        Array<Type> tf(std::move(f));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] -= t;
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator-(const Type &t, Array<Type> &&f) noexcept {
        Array<Type> tf(std::move(f));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = t - tf[i];
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator*(const Array<Type> &f1, const Array<real> &f2) {
        checkArraysSize(f1, f2, "*");
        Array<Type> f(f1.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] * f2[i];
        }

        return f;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator*(const Array<real> &f1, const Array<Type> &f2) {
        checkArraysSize(f1, f2, "*");
        Array<Type> f(f1.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] * f2[i];
        }

        return f;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator*(Array<Type> &&f1, const Array<real> &f2) noexcept {
        checkArraysSize(f1, f2, "*");
        Array<Type> tf(std::move(f1));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] *= f2[i];
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator*(const Array<real> &f1, Array<Type> &&f2) noexcept {
        checkArraysSize(f1, f2, "*");
        Array<Type> tf(std::move(f2));
        for (typename Array<Type>::size_type i = 0; i < f1.size(); ++i) {
            tf[i] = f1[i] * tf[i];
        }

        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator*(const Array<Type> &f, const real &t) {
        Array<Type> tf(f.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            tf[i] = f[i] * t;
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator*(Array<Type> &&f, const real &t) noexcept {
        Array<Type> tf(std::move(f));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] *= t;
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator*(const real t, const Array<Type> &f) {
        Array<Type> tf(f.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            tf[i] = (Type)(t * f[i]);
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard inline Array<Type> operator*(const real t, Array<Type> &&f) noexcept {
        Array<Type> tf(std::move(f));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] = t * tf[i];
        }
        return tf;
    }

    template <class Type>
    hur_nodiscard Array<Type> operator/(const Array<Type> &f1, const Array<real> &f2) {
        checkArraysSize(f1, f2, "/");
        Array<Type> f(f1.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            f[i] = f1[i] / f2[i];
        }

        return f;
    }

    template <class Type>
    hur_nodiscard Array<Type> operator/(Array<Type> &&f1, const Array<real> &f2) noexcept {
        checkArraysSize(f1, f2, "/");
        Array<Type> tf(std::move(f1));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] /= f2[i];
        }
        return tf;
    }

    template <class Type> hur_nodiscard Array<Type> operator/(const Array<Type> &f, const real &t) {
        if (t == real(0)) {
            LFatal("Attempt to devide Array by zero");
        }

        Array<Type> ft(f.size());
        for (typename Array<Type>::size_type i = 0; i < f.size(); ++i) {
            ft[i] = f[i] / t;
        }
        return ft;
    }

    template <class Type>
    hur_nodiscard Array<Type> operator/(Array<Type> &&f, const real &t) noexcept {
        if (t == real(0)) {
            LFatal("Attempt to devide Array by zero");
        }
        Array<Type> tf(std::move(f));
        for (typename Array<Type>::size_type i = 0; i < tf.size(); ++i) {
            tf[i] /= t;
        }
        return tf;
    }

    template <class T1, class T2>
    inline void checkArraysSize(const Array<T1> &l1, const Array<T2> &l2,
                                const std::string operation) {
#ifdef HUR_DEBUG
        if (l1.size() != l2.size()) {
            LFatal("The size of fileds is not equal. Size1 = %d\nSize2 = %d in %s", l1.size(),
                   l2.size(), operation);
        }
#endif // HUR_DEBUG
    }

    template <class Type>
    inline void component(Array<typename Array<Type>::elementType> &comp, const Array<Type> &lf,
                          const int d) {
        checkArraysSize(comp, lf, "component");
        for (typename Array<Type>::size_type i = 0; i < comp.size(); ++i) {
            comp[i] = lf[i].component(d);
        }
    }

    template <class Type> void magSqr(Array<real> &comp, const Array<Type> &lf) {
        checkArraysSize(comp, lf, "magSqr");
        for (typename Array<Type>::size_type i = 0; i < comp.size(); ++i) {
            comp[i] = magSqr(lf[i]);
        }
    }

    template <class Type> hur_nodiscard inline Array<real> magSqr(const Array<Type> &lf) {
        Array<real> MagSqr(lf.size());
        magSqr(MagSqr, lf);
        return MagSqr;
    }

    template <class Type> hur_nodiscard inline Array<real> magSqr(Array<Type> &&lf) {
        Array<real> tf(std::move(lf));
        magSqr(tf, tf);
        return tf;
    }

    template <class Type> void mag(Array<real> &comp, const Array<Type> &lf) {
        checkArraysSize(comp, lf, "mag");
        for (typename Array<Type>::size_type i = 0; i < comp.size(); ++i) {
            comp[i] = mag(lf[i]);
        }
    }

    template <class Type> hur_nodiscard inline Array<real> mag(const Array<Type> &lf) {
        Array<real> Mag(lf.size());
        mag(Mag, lf);
        return Mag;
    }

    template <class Type> hur_nodiscard inline Array<real> mag(Array<Type> &&lf) {
        Array<real> tf(std::move(lf));
        mag(tf, tf);
        return tf;
    }

    template <class Type> hur_nodiscard Type min(const Array<Type> &f) {
        Type minimum = f[0];
        for (auto &bf : f) {
            minimum = min(minimum, bf);
        }
        return minimum;
    }

    template <class Type> hur_nodiscard Type max(const Array<Type> &f) {
        Type maximum = f[0];
        for (auto &bf : f) {
            maximum = max(maximum, bf);
        }
        return maximum;
    }

    template <template <class> class ArrayType, class Type>
    hur_nodiscard inline Type sumArray(ArrayType<Type> &f) {
        Type s = Zero;
        for (typename ArrayType<Type>::size_type i = 0; i < f.size(); ++i) {
            s += f[i];
        }
        return s;
    }
} // namespace OpenHurricane
