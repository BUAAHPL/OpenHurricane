/*!
 * \file FixedList.hpp
 * \brief Header of FixedList
 *       The subroutines and functions are in the <i>FixedList.cpp</i> file.
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

#include "ArrayBase.hpp"
#include "List.hpp"
#include "bools.hpp"
#include "errorAbort.hpp"
#include "zero.hpp"
#include <algorithm>
#include <string>

namespace OpenHurricane {

    // Forward declaration of friend functions and operators

    template <class Type, int nElements> class FixedList;

    template <class Type, int nElements>
    IStringStream &operator>>(IStringStream &, FixedList<Type, nElements> &);

    template <class Type, int nElements>
    OStringStream &operator<<(OStringStream &, const FixedList<Type, nElements> &);

    /**
     * \brief The class of fixed size list.
     */
    template <class Type, int nElements> class FixedList : public ArrayFix<Type, nElements> {
    public:
        using Base = ArrayFix<Type, nElements>;

        using value_type = typename Base::value_type;
        using reference = typename Base::reference;
        using const_reference = typename Base::const_reference;
        using difference_type = typename Base::difference_type;
        using size_type = typename Base::size_type;

        /*!\brief Element Type of value_type.*/
        using elementType = typename feature<value_type>::elementType;

        using iterator = Type *;
        using const_iterator = const Type *;
        using reverse_iterator = Type *;
        using const_reverse_iterator = const Type *;

        static constexpr int nElements_ = Base::nElements_;

    public:
        static const std::streamsize precision;

        hur_nodiscard inline static const FixedList<Type, nElements> &nullObject() {
            return NullRefObj::nullRef<FixedList<Type, nElements>>();
        }

        constexpr inline FixedList() : Base() {}
        constexpr explicit inline FixedList(const Type &t) : Base(t) {}
        constexpr inline FixedList &operator=(const Type &t) {
            Base::operator=(t);
            return *this;
        }

        constexpr explicit inline FixedList(const Type v[nElements]) : Base() {
            for (int i = 0; i < nElements; i++) {
                Base::operator[](i) = v[i];
            }
        }
        constexpr inline FixedList &operator=(const Type v[nElements]) {
            for (int i = 0; i < nElements; i++) {
                Base::operator[](i) = v[i];
            }
            return *this;
        }

        constexpr inline FixedList(std::initializer_list<Type> lst) : Base(lst) {}
        constexpr inline FixedList &operator=(std::initializer_list<Type> lst) {
            Base::operator=(lst);
            return *this;
        }

        constexpr inline FixedList(const FixedList &other) : Base(other) {}
        constexpr inline FixedList &operator=(const FixedList &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        inline FixedList(const List<Type> &other) : Base() {
            for (int i = 0; i < min((int)other.size(), nElements); i++) {
                Base::operator[](i) = other[i];
            }
        }
        inline FixedList &operator=(const List<Type> &other) {
            for (int i = 0; i < min((int)other.size(), nElements); i++) {
                Base::operator[](i) = other[i];
            }
            return *this;
        }

        inline FixedList(IStringStream &is) : Base(is) {}
        inline ~FixedList() noexcept {}

        inline void resize(const size_type) {}

        inline constexpr void transfer(const FixedList &other) {
            for (int i = 0; i < nElements; i++) {
                Base::operator[](i) = other[i];
            }
        }

        hur_nodiscard constexpr inline size_type size() const noexcept { return nElements; }

        hur_nodiscard constexpr inline size_type maxSize() const noexcept { return nElements; }

        hur_nodiscard constexpr inline bool empty() const noexcept { return false; }

        // IOstringstream operators

        friend IStringStream &operator>>
            <Type, nElements>(IStringStream &, FixedList<Type, nElements> &);

        friend OStringStream &operator<< <Type, nElements>(OStringStream &,
                                                           const FixedList<Type, nElements> &);
    };

    template <class Type, int nElements>
    IStringStream &operator>>(IStringStream &is, FixedList<Type, nElements> &L) {
        const std::streamsize defaultPrecision = is.precision();
        is.precision(feature<Type>::precision);
        for (int i = 0; i < nElements; ++i) {
            is >> L[i];
        }

        is.precision(defaultPrecision);
        return is;
    }

    template <class Type, int nElements>
    OStringStream &operator<<(OStringStream &os, const FixedList<Type, nElements> &L) {
        const std::streamsize defaultPrecision = os.precision();
        os.precision(feature<Type>::precision);
        for (int i = 0; i < nElements; ++i) {
            os << L[i] << " ";
        }
        os.precision(defaultPrecision);
        return os;
    }

} // namespace OpenHurricane
