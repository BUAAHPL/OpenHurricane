/*!
 * \file vectorSpace.hpp
 * \brief Header of vector space.
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
#include "ArrayBase.hpp"
#include "HurMPIBase.hpp"
#include "errorAbort.hpp"
#include "logFile.hpp"
#include "real.hpp"
#include "zero.hpp"
#include <string>

namespace OpenHurricane {
    // Forward declaration of friend functions and operators
    template <class Type, class derivedType, int nElements> class VectorSpace;

    template <class Type, class derivedType, int nElements>
    IStringStream &operator>>(IStringStream &, VectorSpace<Type, derivedType, nElements> &);

    template <class Type, class derivedType, int nElements>
    OStringStream &operator<<(OStringStream &, const VectorSpace<Type, derivedType, nElements> &);

    /**
     * \brief The template class of vector space.
     */
    template <class Type, class derivedType, int nElements>
    class VectorSpace : public ArrayFix<Type, nElements> {
    public:
        using Base = ArrayFix<Type, nElements>;

        using vsType = VectorSpace<Type, derivedType, nElements>;

        using value_type = typename Base::value_type;
        using reference = typename Base::reference;
        using const_reference = typename Base::const_reference;
        using difference_type = typename Base::difference_type;
        using size_type = typename Base::size_type;
        using elementType = value_type;

        /**\brief Number of components in this vector space*/
        static constexpr int nElements_ = Base::nElements_;

        /**\brief MPI datatype*/
        static constexpr HurMPIBase::Datatype MPIType = feature<elementType>::MPIType;

        static const std::streamsize precision;

        static constexpr int mRows = nElements;
        static constexpr int nCols = 1;

    public:

        constexpr inline VectorSpace() : Base() {}

        constexpr inline VectorSpace(const OpenHurricane::zero) : Base(Zero) {}
        constexpr inline VectorSpace &operator=(const OpenHurricane::zero) {
            Base::operator=(Zero);
            return *this;
        }

        constexpr inline VectorSpace(const VectorSpace &other) : Base(other) {}
        constexpr inline VectorSpace &operator=(const VectorSpace &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        template <class otherType, class otherDerived>
        constexpr inline explicit VectorSpace(
            const VectorSpace<otherType, otherDerived, nElements> &other)
            : Base(other) {}
        template <class otherType, class otherDerived>
        constexpr inline vsType &
        operator=(const VectorSpace<otherType, otherDerived, nElements> &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        constexpr inline VectorSpace(const value_type &elem) : Base(elem) {}

        inline VectorSpace(IStringStream &is) : Base() {
            std::string str;
            is >> str;
            replaceAllMarks(str, "(", " ");
            replaceAllMarks(str, ")", " ");
            replaceAllMarks(str, ",", " ");
            IStringStream sstr(str);
            sstr.precision(feature<value_type>::precision);
            for (integer i = 0; i < nElements; ++i) {
                sstr >> Base::operator[](i);
            }
        }

        inline constexpr VectorSpace(std::initializer_list<value_type> lst) : Base(lst) {}
        inline constexpr VectorSpace &operator=(std::initializer_list<value_type> lst) {
            Base::operator=(lst);
            return *this;
        }

        inline ~VectorSpace() noexcept {}

        hur_nodiscard constexpr inline const value_type &
        component(const size_type i) const noexcept {
            return Base::operator[](i);
        }

        hur_nodiscard constexpr inline value_type &component(const size_type i) noexcept {
            return Base::operator[](i);
        }

        constexpr inline void component(value_type &elem, const size_type i) const noexcept {
            Base::element(elem, i);
        }

        constexpr inline void setConstant(const value_type &elem) { Base::setConstant(elem); }
        constexpr inline void setZero() { Base::setZero(); }

        /**\brief Return the magnitude of VectorSpace*/
        hur_nodiscard inline value_type magnitude() const {
            value_type mag = value_type(0);
            for (int i = 0; i < nElements_; i++) {
                mag += Base::operator[](i) * Base::operator[](i);
            }
            return value_type(sqrt(mag));
        }

        /**\brief Return the magnitude square of VectorSpace*/
        hur_nodiscard constexpr inline value_type magSqr() const {
            value_type mag = value_type(0);
            for (int i = 0; i < nElements_; i++) {
                mag += Base::operator[](i) * Base::operator[](i);
            }
            return mag;
        }

        /**\brief Return a normalized vector of VectorSpace*/
        hur_nodiscard inline VectorSpace normalized() const {
            if (this->operator==(OpenHurricane::Zero)) {
                LFatal("Attempt to normalize a zero vectorspace!");
            }
            vsType vt;
            value_type mag = this->magnitude();
            for (int i = 0; i < nElements_; i++) {
                vt[i] = Base::operator[](i) / mag;
            }
            return vt;
        }

        constexpr inline void operator+=(const VectorSpace &other) noexcept {
            for (int i = 0; i < nElements_; i++) {
                Base::operator[](i) += other[i];
            }
        }

        constexpr inline void operator-=(const VectorSpace &other) noexcept {
            for (int i = 0; i < nElements_; i++) {
                Base::operator[](i) -= other[i];
            }
        }

        hur_nodiscard constexpr inline VectorSpace
        operator+(const VectorSpace &other) const noexcept {
            vsType v2;
            for (int i = 0; i < nElements_; i++) {
                v2[i] = Base::operator[](i) + other[i];
            }
            return v2;
        }

        hur_nodiscard constexpr inline VectorSpace
        operator-(const VectorSpace &other) const noexcept {
            vsType v2;
            for (int i = 0; i < nElements_; i++) {
                v2[i] = Base::operator[](i) - other[i];
            }
            return v2;
        }

        constexpr inline void operator*=(const value_type &elem) noexcept {
            for (int i = 0; i < nElements_; i++) {
                Base::operator[](i) *= elem;
            }
        }

        constexpr inline void operator/=(const value_type &elem) {
            if (elem == value_type(0)) {
                LFatal("divided by zero ! elem = %s", std::to_string(elem).c_str());
            }
            for (int i = 0; i < nElements_; i++) {
                Base::operator[](i) /= elem;
            }
        }

        hur_nodiscard constexpr inline bool operator==(const VectorSpace &other) const noexcept {
            bool equal = true;
            for (int i = 0; i < nElements_; i++) {
                equal = (equal && (Base::operator[](i) == other[i]));
            }
            return equal;
        }

        hur_nodiscard constexpr inline bool operator!=(const VectorSpace &other) const noexcept {
            return !this->operator==(other);
        }

        hur_nodiscard constexpr inline bool operator==(const OpenHurricane::zero) const noexcept {
            bool equal = true;
            for (int i = 0; i < nElements_; i++) {
                equal = (equal && (Base::operator[](i) == value_type(OpenHurricane::Zero)));
            }
            return equal;
        }

        hur_nodiscard constexpr inline bool operator!=(const OpenHurricane::zero) const noexcept {
            return !this->operator==(OpenHurricane::Zero);
        }

        friend IStringStream &operator>>
            <Type, derivedType, nElements>(IStringStream &,
                                           VectorSpace<Type, derivedType, nElements> &);

        friend OStringStream &operator<< <Type, derivedType, nElements>(
            OStringStream &, const VectorSpace<Type, derivedType, nElements> &);
    };

    template <class Type, class derivedType, int nElements>
    IStringStream &operator>>(IStringStream &is, VectorSpace<Type, derivedType, nElements> &L) {
        std::string str;
        is >> str;
        replaceAllMarks(str, "(", " ");
        replaceAllMarks(str, ")", " ");
        replaceAllMarks(str, ",", " ");
        IStringStream sstr(str);
        sstr.precision(feature<Type>::precision);
        for (integer i = 0; i < nElements; ++i) {
            sstr >> L[i];
        }
        return is;
    }

    template <class Type, class derivedType, int nElements>
    OStringStream &operator<<(OStringStream &os,
                              const VectorSpace<Type, derivedType, nElements> &L) {
        const std::streamsize defaultPrecision = os.precision();
        os.precision(feature<Type>::precision);
        os << "(";
        for (integer i = 0; i < nElements; ++i) {
            os << L[i];
            if (i < nElements - 1) {
                os << ",";
            }
        }
        os << ")";
        return os;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    operator/(const VectorSpace<Type, derivedType, nElements> &vs, real s) {
        derivedType v;
        for (int i = 0; i < nElements; i++) {
            v[i] = vs[i] / s;
        }
        return v;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentAdd(const VectorSpace<Type, derivedType, nElements> &vs1, const Type &v2) {
        derivedType v;
        for (int i = 0; i < nElements; ++i) {
            v[i] = componentAdd(vs1[i], v2);
        }
        return v;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentAdd(const Type &v1, const VectorSpace<Type, derivedType, nElements> &vs2) {
        derivedType v;
        for (int i = 0; i < nElements; ++i) {
            v[i] = componentAdd(v1, vs2[i]);
        }
        return v;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentSubtract(const VectorSpace<Type, derivedType, nElements> &vs1, const Type &v2) {
        derivedType v;

        for (int i = 0; i < nElements; ++i) {
            v[i] = componentSubtract(vs1[i], v2);
        }
        return v;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentSubtract(const Type &v1, const VectorSpace<Type, derivedType, nElements> &vs2) {
        derivedType v;

        for (int i = 0; i < nElements; ++i) {
            v[i] = componentSubtract(v1, vs2[i]);
        }
        return v;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentMultiply(const VectorSpace<Type, derivedType, nElements> &vs1,
                      const VectorSpace<Type, derivedType, nElements> &vs2) {
        derivedType v;

        for (int i = 0; i < nElements; ++i) {
            v[i] = componentMultiply(vs1[i], vs2[i]);
        }
        return v;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard hur_forceinline derivedType
    componentPow(const VectorSpace<Type, derivedType, nElements> &vs1,
                 const VectorSpace<Type, derivedType, nElements> &vs2) {
        derivedType v;

        for (int i = 0; i < nElements; ++i) {
            v[i] = componentPow(vs1[i], vs2[i]);
        }
        return v;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentDivide(const VectorSpace<Type, derivedType, nElements> &vs1,
                    const VectorSpace<Type, derivedType, nElements> &vs2) {
        derivedType v;

        for (int i = 0; i < nElements; ++i) {
            v[i] = componentDivide(vs1[i], vs2[i]);
        }
        return v;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline Type
    componentMax(const VectorSpace<Type, derivedType, nElements> &vs) {
        Type s = vs[0];

        for (int i = 1; i < nElements; ++i) {
            s = max(s, vs[i]);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentMax(const VectorSpace<Type, derivedType, nElements> &vs1,
                 const VectorSpace<Type, derivedType, nElements> &vs2) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = max(vs1[i], vs2[i]);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentMax(const Type s1, const VectorSpace<Type, derivedType, nElements> &vs2) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = max(s1, vs2[i]);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentMax(const VectorSpace<Type, derivedType, nElements> &vs1, const Type s2) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = max(vs1[i], s2);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline Type
    componentMin(const VectorSpace<Type, derivedType, nElements> &vs) {
        Type s = vs[0];

        for (int i = 1; i < nElements; ++i) {
            s = min(s, vs[i]);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentMin(const VectorSpace<Type, derivedType, nElements> &vs1,
                 const VectorSpace<Type, derivedType, nElements> &vs2) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = min(vs1[i], vs2[i]);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentMin(const Type s1, const VectorSpace<Type, derivedType, nElements> &vs2) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = min(s1, vs2[i]);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentMin(const VectorSpace<Type, derivedType, nElements> &vs1, const Type s2) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = min(vs1[i], s2);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline Type
    componentAv(const VectorSpace<Type, derivedType, nElements> &vs) {
        Type s = vs[0];

        for (int i = 1; i < nElements; ++i) {
            s += vs[i];
        }
        return s / nElements;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentSqr(const VectorSpace<Type, derivedType, nElements> &vs) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = sqr(vs[i]);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard inline derivedType
    componentSqrt(const VectorSpace<Type, derivedType, nElements> &vs) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = sqrt(vs[i]);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard inline derivedType
    componentMag(const VectorSpace<Type, derivedType, nElements> &vs) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = mag(vs[i]);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    componentSign(const VectorSpace<Type, derivedType, nElements> &vs) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = componentSign(vs[i]);
        }
        return s;
    }

    template <class Type, class derivedType, int nElements>
    hur_nodiscard constexpr inline derivedType
    inv(const VectorSpace<Type, derivedType, nElements> &vs) {
        derivedType s;

        for (int i = 0; i < nElements; ++i) {
            s[i] = inv(vs[i]);
        }
        return s;
    }

} // namespace OpenHurricane

#include "productTypes.hpp"
