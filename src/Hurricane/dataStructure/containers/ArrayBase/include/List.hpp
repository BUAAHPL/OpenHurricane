/*!
 * \file List.hpp
 * \brief Header of List.
 *       The subroutines and functions are in the <i>List.inl</i> file.
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
#include "basicFunctions.hpp"
#include "bools.hpp"
#include "errorAbort.hpp"
#include "integer.hpp"
#include "smartPointer.hpp"
#include "zero.hpp"
#include <algorithm>
#include <initializer_list>
#include <random>
#include <stdarg.h>
#include <string>

namespace OpenHurricane {
    // Forward declaration of friend functions and operators

    template <class Type> class List;

    template <class Type> IStringStream &operator>>(IStringStream &, List<Type> &);

    template <class Type> OStringStream &operator<<(OStringStream &, const List<Type> &);

    /**
     * \brief The class of list.
     */
    template <class Type> class List : public ArrayDyn<Type, integer> {
    public:
        using Base = ArrayDyn<Type, integer>;

        using value_type = typename Base::value_type;
        using reference = typename Base::reference;
        using const_reference = typename Base::const_reference;
        using difference_type = typename Base::difference_type;
        using size_type = typename Base::size_type;

    public:
        static constexpr size_type npos = Base::npos;

        /**\brief Return a null List.*/
        hur_nodiscard inline static const List<Type> &nullObject() {
            return NullRefObj::nullRef<List<Type>>();
        }

        inline List() : Base() {}
        inline explicit List(const size_type s) : Base(s) {}

        inline List(const size_type s, const Type &val) : Base(s, val) {}
        inline List(const size_type s, const Type &val, const size_type capa)
            : Base(s, val, capa) {}
        inline List &operator=(const Type &elem) {
            Base::operator=(elem);
            return *this;
        }

        /**\brief Construct with given size initializing all elements to zero.*/
        inline List(const size_type s, const zero) : Base(s, Zero) {}
        inline List &operator=(const zero) {
            Base::operator=(Zero);
            return *this;
        }

        /**\brief Copy constructor.*/
        inline List(const List &other) : Base(other) {}
        inline List &operator=(const List &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }
        inline List(const Base &other) : Base(other) {}
        inline List &operator=(const Base &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        inline List(List &&other) noexcept : Base(std::move(other)) {}
        inline List &operator=(List &&other) noexcept {
            Base::operator=(std::move(other));
            return *this;
        }

        template <class otherType>
        inline explicit List(const List<otherType> &other) : Base(other) {}
        template <class otherType> List &operator=(const List<otherType> &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        inline List(std::initializer_list<Type> lst) : Base(lst) {}
        inline List &operator=(std::initializer_list<Type> lst) {
            Base::operator=(lst);
            return *this;
        }

        inline List(IStringStream &is) : Base(is) {}

        hur_nodiscard inline uniquePtr<List> clone() const {
            return uniquePtr<List>(new List(*this));
        }

        /**\brief Destructor.*/
        inline virtual ~List() noexcept {}

        /*!\brief Return a sub-array.*/
        hur_nodiscard inline List sub(const size_type Off = 0, const size_type Count = npos) {
            return List(Base::sub(Off, Count));
        }

        // IOstringstream operators

        friend IStringStream &operator>><Type>(IStringStream &, List &);

        friend OStringStream &operator<< <Type>(OStringStream &, const List &);
    };

    template <class Type> void sort(List<Type> &a) {
        std::sort(a.begin(), a.end());
    }

    template <class Type, class Cmp> void sort(List<Type> &a, const Cmp &cmp) {
        std::sort(a.begin(), a.end(), cmp);
    }

    template <class Type> void stableSort(List<Type> &a) {
        std::stable_sort(a.begin(), a.end());
    }

    template <class Type, class Cmp> void stableSort(List<Type> &a, const Cmp &cmp) {
        std::stable_sort(a.begin(), a.end(), cmp);
    }

    template <class Type> void shuffle(List<Type> &a) {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(a.begin(), a.end(), g);
    }

    // Reverse the first n elements of the list
    template <class Type>
    inline void reverse(List<Type> &a, const typename List<Type>::size_type n) {
        for (typename List<Type>::size_type i = 0; i < n / 2; i++) {
            Swap(a[i], a[n - 1 - i]);
        }
    }

    // Reverse all the elements of the list
    template <class Type> inline void reverse(List<Type> &a) {
        reverse(a, a.size());
    }

    template <class Type> inline std::string toString(const List<Type> &l) {
        std::string str;
        str = "";
        for (typename List<Type>::size_type i = 0; i < l.size(); i++) {
            str += OpenHurricane::toString(l[i]);
            if ((i > 0) && (i % 10 == 0)) {
                str.append("\n");
                continue;
            }
            str += " ";
        }
        return str;
    }

    template <class Type> IStringStream &operator>>(IStringStream &is, List<Type> &L) {
        const std::streamsize defaultPrecision = is.precision();
        typename List<Type>::size_type Size;

        is >> Size;
        L.resize(Size);

        is.precision(feature<Type>::precision);
        for (typename List<Type>::size_type i = 0; i < Size; ++i) {
            is >> L[i];
        }

        is.precision(defaultPrecision);
        return is;
    }

    template <class Type> OStringStream &operator<<(OStringStream &os, const List<Type> &L) {
        const std::streamsize defaultPrecision = os.precision();

        const typename List<Type>::size_type Size = L.size();
        os << Size << " ";

        os.precision(feature<Type>::precision);
        for (typename List<Type>::size_type i = 0; i < Size; ++i) {
            os << L[i] << " ";
        }
        os.precision(defaultPrecision);
        return os;
    }

} // namespace OpenHurricane
