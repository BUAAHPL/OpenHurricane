/*!
 * \file ArrayBase.hpp
 * \brief Headers of the base of Array.
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

#include "ArrayStorage.hpp"
#include "NullRefObj.hpp"

namespace OpenHurricane {

    template <class Type, int nElements> class ArrayFix {
    public:
        using Base = storageFix<Type, nElements>;

        using value_type = typename Base::value_type;
        using reference = typename Base::reference;
        using const_reference = typename Base::const_reference;
        using difference_type = typename Base::difference_type;
        using size_type = typename Base::size_type;

        using iterator = value_type *;
        using const_iterator = const value_type *;
        using reverse_iterator = value_type *;
        using const_reverse_iterator = const value_type *;

        static constexpr int nElements_ = Base::nElements_;

    private:
        Base mData_;

    public:
        hur_nodiscard inline static const ArrayFix &nullObject() {
            return NullRefObj::nullRef<ArrayFix>();
        }

        static constexpr size_type npos = -2;

        /**
         * \brief Null constructor.
         */
        constexpr inline ArrayFix() : mData_() {}

        /**
         * \brief Construct with given size and value for all elements.
         */
        constexpr inline ArrayFix(const value_type &elem) : mData_(elem) {}
        constexpr inline ArrayFix &operator=(const value_type &elem) {
            for (size_type i = 0; i < this->size(); ++i) {
                *(begin() + i) = elem;
            }
            return *this;
        }

        /**
         * \brief Construct with given size and zero for all elements.
         */
        constexpr inline ArrayFix(const zero) : mData_(Zero) {}
        constexpr inline ArrayFix &operator=(const zero) {
            for (size_type i = 0; i < this->size(); ++i) {
                *(begin() + i) = Zero;
            }
            return *this;
        }

        /**
         * \brief Copy constructor.
         */
        constexpr inline ArrayFix(const ArrayFix &other) : mData_(other.mData_) {}
        constexpr inline ArrayFix &operator=(const ArrayFix &other) {
            if (this != std::addressof(other)) {
                mData_ = other.mData_;
            }
            return *this;
        }

        template <class otherType>
        constexpr inline ArrayFix(const ArrayFix<otherType, nElements> &other)
            : mData_(other.mData_) {}
        template <class otherType>
        constexpr inline ArrayFix &operator=(const ArrayFix<otherType, nElements> &other) {
            if (this != std::addressof(other)) {
                mData_ = other.mData_;
            }
            return *this;
        }

        constexpr ArrayFix(std::initializer_list<value_type> lst) : mData_() {
            if (mData_.size() > 0) {
                auto iter = lst.begin();
                for (size_type i = 0; i < min(mData_.size(), (integer)lst.size()); ++i) {
                    *(mData_.data() + i) = *(iter++);
                }
            }
        }
        constexpr ArrayFix &operator=(std::initializer_list<value_type> lst) {
            size_type sizes = size_type(std::distance(lst.begin(), lst.end()) + 1);

            if (mData_.size() != 0) {
                auto iter = lst.begin();
                for (size_type i = 0; i < min(mData_.size(), (integer)lst.size()); ++i) {
                    *(mData_.data() + i) = *(iter++);
                }
            }
            return *this;
        }

        ArrayFix(IStringStream &is) : mData_() {
            const std::streamsize defaultPrecision = is.precision();
            is.precision(feature<value_type>::precision);
            size_type Size = 0;
            is >> Size;
            for (size_type i = 0; i < min(mData_.size(), Size); ++i) {
                is >> *(mData_.data() + i);
            }
            is.precision(defaultPrecision);
        }

        virtual inline ~ArrayFix() noexcept {}

        inline void swap(ArrayFix &other) noexcept { mData_.swap(other.mData_); }

        constexpr inline void setConstant(const value_type &elem) {
            for (size_type i = 0; i < this->size(); ++i) {
                *(begin() + i) = elem;
            }
        }

        constexpr inline void setZero() {
            for (size_type i = 0; i < this->size(); ++i) {
                *(begin() + i) = Zero;
            }
        }

        /**
         * \brief Return the number of elements in the Array.
         */
        hur_nodiscard inline constexpr size_type size() const noexcept { return mData_.size(); }

        hur_nodiscard inline constexpr const value_type *cdata() const noexcept {
            return mData_.data();
        }
        hur_nodiscard inline constexpr const value_type *data() const noexcept {
            return mData_.data();
        }
        hur_nodiscard constexpr inline value_type *data() noexcept { return mData_.data(); }

        hur_nodiscard constexpr inline value_type &first() noexcept { return mData_[0]; }

        hur_nodiscard inline constexpr const value_type &first() const noexcept {
            return data()[0];
        }
        hur_nodiscard inline constexpr value_type &front() noexcept { return data()[0]; }
        hur_nodiscard inline constexpr const value_type &front() const noexcept {
            return data()[0];
        }
        hur_nodiscard constexpr inline value_type &last() noexcept { return data()[size() - 1]; }
        hur_nodiscard inline constexpr const value_type &last() const noexcept {
            return data()[size() - 1];
        }
        hur_nodiscard constexpr inline value_type &back() noexcept { return data()[size() - 1]; }
        hur_nodiscard inline constexpr const value_type &back() const noexcept {
            return data()[size() - 1];
        }

        hur_nodiscard constexpr inline iterator begin() noexcept { return mData_.data(); }
        hur_nodiscard constexpr inline iterator end() noexcept { return (mData_.data() + size()); }

        hur_nodiscard constexpr inline const_iterator cbegin() const noexcept {
            return mData_.data();
        }
        hur_nodiscard constexpr inline const_iterator cend() const noexcept {
            return (mData_.data() + size());
        }

        hur_nodiscard constexpr inline const_iterator begin() const noexcept {
            return mData_.data();
        }
        hur_nodiscard constexpr inline const_iterator end() const noexcept {
            return (mData_.data() + size());
        }

        hur_nodiscard constexpr inline reverse_iterator rbegin() noexcept {
            return (mData_.data() + size() - 1);
        }
        hur_nodiscard constexpr inline reverse_iterator rend() noexcept {
            return mData_.data() - 1;
        }

        hur_nodiscard constexpr inline const_reverse_iterator crbegin() const noexcept {
            return (mData_.data() + size() - 1);
        }
        hur_nodiscard constexpr inline const_reverse_iterator crend() const noexcept {
            return mData_.data() - 1;
        }

        hur_nodiscard constexpr inline const_reverse_iterator rbegin() const noexcept {
            return (mData_.data() + size() - 1);
        }
        hur_nodiscard constexpr inline const_reverse_iterator rend() const noexcept {
            return mData_.data() - 1;
        }

        constexpr inline void element(value_type &ele, const size_type i) const noexcept {
            mData_.element(ele, i);
        }

        constexpr inline void replace(const size_type i, const value_type &ele) noexcept {
            mData_.replace(i, ele);
        }

        hur_nodiscard constexpr inline value_type &operator[](const size_type i) noexcept {
            return mData_.element(i);
        }

        hur_nodiscard constexpr inline const value_type &
        operator[](const size_type i) const noexcept {
            return mData_.element(i);
        }

        /**\brief Return element of List.*/
        hur_nodiscard constexpr inline value_type &operator()(const size_type i) noexcept {
            return mData_.element(i);
        }

        /**
         *\brief Return element of constant List
         */
        hur_nodiscard constexpr inline const value_type &
        operator()(const size_type i) const noexcept {
            return mData_.element(i);
        }

        hur_nodiscard constexpr bool operator==(const ArrayFix &other) const noexcept {
            for (size_type i = 0; i < this->size(); ++i) {
                if (*(begin() + i) != *(other.begin() + i)) {
                    return false;
                }
            }
            return true;
        }

        hur_nodiscard constexpr inline bool operator!=(const ArrayFix &other) const noexcept {
            return !operator==(other);
        }
    };

    template <class Type, int nElements> void sort(ArrayFix<Type, nElements> &arr) {
        std::sort(arr.begin(), arr.end());
    }
    template <class Type, int nElements, class compare>
    void sort(ArrayFix<Type, nElements> &arr, const compare &cmp) {
        std::sort(arr.begin(), arr.end(), cmp);
    }

    template <class Type, int nElements> void stableSort(ArrayFix<Type, nElements> &arr) {
        std::stable_sort(arr.begin(), arr.end());
    }

    template <class Type, int nElements, class compare>
    void stableSort(ArrayFix<Type, nElements> &ar, const compare &cmp) {
        std::stable_sort(ar.begin(), ar.end(), cmp);
    }

    template <class Type, int nElements>
    hur_nodiscard inline std::string toString(const ArrayFix<Type, nElements> &arr) {
        std::string str;
        str = "";
        integer i = 0;
        for (const auto &e : arr) {
            str += toString(e);
            if (i < arr.size() - 1) {
                str += ", ";
            }
            i++;
            if ((i % 10 == 0)) {
                str.append("\n");
            }
        }
        return str;
    }

    template <class Type, typename sizeType> class ArrayDyn {
    public:
        using Base = storageDyn<Type, sizeType>;

        using value_type = typename Base::value_type;
        using reference = typename Base::reference;
        using const_reference = typename Base::const_reference;
        using difference_type = typename Base::difference_type;
        using size_type = typename Base::size_type;

        using iterator = value_type *;
        using const_iterator = const value_type *;
        using reverse_iterator = value_type *;
        using const_reverse_iterator = const value_type *;

    private:
        Base mData_;

    public:
        hur_nodiscard inline static const ArrayDyn &nullObject() {
            return NullRefObj::nullRef<ArrayDyn>();
        }

        static constexpr size_type npos = -2;

        /**
         * \brief Null constructor.
         */
        inline ArrayDyn() : mData_() {}

        /**
         * \brief Construct with given size.
         */
        inline explicit ArrayDyn(const size_type newSize) : mData_(newSize) {}

        /**
         * \brief Construct with given size and value for all elements.
         */
        inline ArrayDyn(const size_type newSize, const value_type &elem)
            : mData_(newSize, newSize, elem) {}
        inline ArrayDyn &operator=(const value_type &elem) {
            for (size_type i = 0; i < this->size(); ++i) {
                *(begin() + i) = elem;
            }
            return *this;
        }

        /**
         * \brief Construct with given size and zero for all elements.
         */
        inline ArrayDyn(const size_type newSize, const zero) : mData_(newSize, newSize, Zero) {}
        inline ArrayDyn &operator=(const zero) {
            for (size_type i = 0; i < this->size(); ++i) {
                *(begin() + i) = Zero;
            }
            return *this;
        }

        /**
         * \brief Construct with given size, capacity and value for all elements.
         */
        inline ArrayDyn(const size_type newSize, const value_type &elem, const size_type newCpacity)
            : mData_(newSize, newCpacity, elem) {}

        /**
         * \brief Construct with given size, capacity and zero for all elements.
         */
        inline ArrayDyn(const size_type newSize, const zero, const size_type newCpacity)
            : mData_(newSize, newCpacity, Zero) {}

        /**
         * \brief Copy constructor.
         */
        inline ArrayDyn(const ArrayDyn &other) : mData_(other.mData_) {}
        inline ArrayDyn &operator=(const ArrayDyn &other) {
            if (this != std::addressof(other)) {
                mData_ = other.mData_;
            }
            return *this;
        }

        /**
         * \brief Construct as copy from rvalue.
         */
        inline ArrayDyn(ArrayDyn &&other) noexcept : mData_(std::move(other.mData_)) {}
        inline ArrayDyn &operator=(ArrayDyn &&other) {
            if (this != std::addressof(other)) {
                mData_.swap(other.mData_);
            }
            return *this;
        }

        template <class otherType>
        inline ArrayDyn(const ArrayDyn<otherType, sizeType> &other) : mData_(other.size()) {
            for (size_type i = 0; i < this->size(); ++i) {
                *(mData_.data() + i) = other[i];
            }
        }
        template <class otherType> ArrayDyn &operator=(const ArrayDyn<otherType, sizeType> &other) {
            if (this != std::addressof(other)) {
                mData_.resize(other.size());
                for (size_type i = 0; i < this->size(); ++i) {
                    *(mData_.data() + i) = other[i];
                }
            }
            return *this;
        }

        ArrayDyn(std::initializer_list<value_type> lst)
            : mData_(static_cast<size_type>(lst.size())) {
            if (mData_.size() > 0) {
                auto iter = lst.begin();
                for (size_type i = 0; i < mData_.size(); ++i) {
                    *(mData_.data() + i) = *(iter++);
                }
            }
        }
        ArrayDyn &operator=(std::initializer_list<value_type> lst) {
            size_type sizes = size_type(std::distance(lst.begin(), lst.end()) + 1);
            resize(sizes);
            if (mData_.size() != 0) {
                auto iter = lst.begin();
                for (size_type i = 0; i < mData_.size(); ++i) {
                    *(mData_.data() + i) = *(iter++);
                }
            }
            return *this;
        }

        ArrayDyn(IStringStream &is) : mData_() {
            const std::streamsize defaultPrecision = is.precision();
            is.precision(feature<value_type>::precision);
            size_type Size = 0;
            is >> Size;
            mData_.resize(Size);
            for (size_type i = 0; i < mData_.size(); ++i) {
                is >> *(mData_.data() + i);
            }
            is.precision(defaultPrecision);
        }

        virtual inline ~ArrayDyn() noexcept {}

        /*!\brief Return a sub-array.*/
        hur_nodiscard ArrayDyn<Type, sizeType> sub(const size_type Off = 0,
                                                   const size_type Count = npos) {
            if (Off < 0 || (Off != 0 && Off >= mData_.size())) {
                LFatal("Offset: %d is out of range 0 ... %d", Off, max(mData_.size() - 1, 0));
            }
            size_type subSize;
            if (Count == npos) {
                subSize = mData_.size() - Off;
            } else {
                if (Count < 0 || (Count != 0 && Count > mData_.size())) {
                    LFatal("Count: %d is out of range 0 ... %d", Count, max(mData_.size() - 1, 0));
                }
                if (Count > (mData_.size() - Off)) {
                    subSize = mData_.size() - Off;
                } else {
                    subSize = Count;
                }
            }
            ArrayDyn<Type, sizeType> subArr(subSize);
            for (size_type i = 0; i < subArr.size(); i++) {
                subArr[i] = *(mData_.data() + i + Off);
            }
            return subArr;
        }

        /**
         * \brief Clear the Array, i.e. set size to zero..
         */
        inline void clear() noexcept { mData_.clear(); }

        inline void swap(ArrayDyn &other) noexcept {
            if (this != std::addressof(other)) {
                mData_.swap(other.mData_);
            }
        }

        /**
         * \brief Transfer the elements of the argument Array into this Array and clear the argument Array.
         */
        inline void transfer(ArrayDyn &other) noexcept {
            this->swap(other);
            other.clear();
        }

        /**
         * \brief Reset size of Array. Only for dynamic Array.
         */
        inline void resize(const size_type newSize) { mData_.resize(newSize); }

        /**
         * \brief Reset size of Array and value for new elements. Only for dynamic Array.
         */
        inline void resize(const size_type newSize, const value_type &elem) {
            mData_.resize(newSize, elem);
        }

        inline void resize(const size_type newSize, const zero) {
            mData_.resize(newSize);
            setZero();
        }

        /**
         * \brief Enlarge the capacity of Array. Only for dynamic Array.
         */
        inline void reserve(const size_type newCapacity) { mData_.reserve(newCapacity); }

        /**
         * \brief Append an element at the end of the Array. Only for dynamic Array.
         */
        inline void append(const value_type &elem) {
            mData_.resize(size() + 1);
            this->last() = elem;
        }

        /**
         * \brief Append a Array at the end of this Array. Only for dynamic Array.
         */
        inline void append(const ArrayDyn<value_type, sizeType> &lst) {
            auto oldSize = mData_.size();
            size_type newSize = oldSize + lst.size();
            mData_.resize(newSize);
            if (mData_.size() != oldSize) {
                for (size_type i = 0; i < lst.size(); i++) {
                    mData_[oldSize++] = lst[i];
                }
            }
        }

        /**
         * \brief Append a Array at the end of this Array. Only for dynamic Array.
         */
        template <int sizeOfElement>
        inline void append(const ArrayFix<value_type, sizeOfElement> &lst) {
            auto oldSize = mData_.size();
            size_type newSize = oldSize + lst.size();
            mData_.resize(newSize);
            if (mData_.size() != oldSize) {
                for (size_type i = 0; i < lst.size(); i++) {
                    mData_[oldSize++] = lst[i];
                }
            }
        }

        /**
         * \brief Insert an element at the entry specified by iter of the Array. Only for dynamic Array.
         * \param[in] iter - The entry where to insert.
         * \param[in] elem - The value to be inserted.
         * \return Iterator pointing to the inserted elem.
         */
        inline iterator insert(const_iterator iter, const value_type &elem) {
            return mData_.insert(iter, elem);
        }

        /**
         * \brief Insert an element at the entry specified by iter of the Array. Only for dynamic Array.
         * \param[in] iter - The entry where to insert.
         * \param[in] num - Specify the total number to be insert.
         * \param[in] elem - The value to be inserted.
         * \return Iterator pointing to the first element inserted.
         */
        inline iterator insert(const_iterator iter, const size_type num, const value_type &elem) {
            return mData_.insert(iter, num, elem);
        }

        /**
         * \brief Insert an Array at the entry specified by iter of the Array. Only for dynamic Array.
         * \param[in] iter - The entry where to insert.
         * \param[in] arr - The Array to be insert.
         * \return Iterator pointing to the first element inserted.
         */
        inline iterator insert(const_iterator iter, const ArrayDyn<value_type, sizeType> &arr) {
            if ((iter > this->end()) || (iter < this->begin())) {
                LFatal("Beyond the size, please check!");
            }
            if (arr.size() > 0) {
                size_type newSize = this->size() + arr.size();
                Base tmpArray(newSize, max(this->capacity(), newSize));

                value_type *oldElem = mData_.data();
                value_type *newElem = tmpArray.data();
                size_type i = this->size();
                iterator newIter = begin();
                while (i--) {
                    if (oldElem == iter) {
                        for (size_type j = 0; j < arr.size(); j++) {
                            if (j == 0) {
                                newIter = newElem;
                            }
                            *newElem++ = arr[j];
                        }
                    }
                    *newElem++ = *oldElem++;
                }
                if (iter == this->end()) {
                    for (size_type j = 0; j < arr.size(); j++) {
                        if (j == 0) {
                            newIter = newElem;
                        }
                        *newElem++ = arr[j];
                    }
                }
                mData_.swap(tmpArray);
                return newIter;
            }
            return const_cast<iterator>(iter);
        }

        /**
         * \brief Erase an element at the entry specified by iter of the Array. Only for dynamic Array.
         * \param[in] iter - The entry where to erase.
         */
        inline iterator erase(const_iterator iter) {
            if ((iter > this->end()) || (iter < this->begin())) {
                LFatal("Beyond the size, please check!");
            }
            size_type newSize = this->size() - 1;
            Base tmpArray(newSize, this->capacity());

            value_type *oldElem = mData_.data();
            value_type *newElem = tmpArray.data();
            size_type i = this->size();
            bool erased = false;
            bool isEnd = false;
            iterator nextIter = begin();
            while (i--) {
                if (oldElem == iter && !erased) {
                    oldElem++;
                    i--;
                    erased = true;
                }
                if (i < 0) {
                    if (erased) {
                        isEnd = true;
                    }
                    break;
                }
                if (erased) {
                    nextIter = newElem;
                    erased = false;
                }
                *newElem++ = *oldElem++;
            }
            mData_.swap(tmpArray);
            if (isEnd) {
                nextIter = this->end();
            }
            return nextIter;
        }

        /**
         * \brief Erase num elements at the entry specified by iter of the Array. Only for dynamic Array.
         * \param[in] iter - The entry where to erase.
         * \param[in] num - The number of elements to be erase.
         * \return Iterator following the last removed element.
         */
        inline iterator erase(const_iterator iter, const size_type num) {
            if ((iter > this->end()) || (iter < this->begin())) {
                LFatal("Beyond the size, please check!");
            }
            size_type newSize = this->size() - 1;
            Base tmpArray(newSize, this->capacity());

            value_type *oldElem = mData_.data();
            value_type *newElem = tmpArray.data();
            size_type i = this->size();
            bool erased = false;
            bool isEnd = false;
            iterator nextIter = begin();
            while (i--) {
                if (oldElem == iter && !erased) {
                    for (size_type j = 0; j < num; j++) {
                        oldElem++;
                        i--;
                    }
                }
                if (i < 0) {
                    if (erased) {
                        isEnd = true;
                    }
                    break;
                }
                if (erased) {
                    nextIter = newElem;
                    erased = false;
                }
                *newElem++ = *oldElem++;
            }
            mData_.swap(tmpArray);
            if (isEnd) {
                nextIter = this->end();
            }
            return nextIter;
        }

        /**
         * \brief Erase the last element of the Array.
         */
        inline void pop_back() { erase(&mData_[this->size() - 1]); }

        /**
         * \brief Erase the first element of the Array.
         */
        inline void pop_front() { erase(this->begin()); }

        /**
         * \brief Insert an element at the last entry of the Array.
         */
        inline void push_back(const value_type &elem) { insert(this->end(), elem); }

        /**
         * \brief Insert an element at the begin entry of the Array.
         */
        inline void push_front(const value_type &elem) { insert(this->begin(), elem); }

        /**
         * \brief Remove all the entries that are same with value of the Array.
         */
        inline void remove(const value_type &elem) { mData_.remove(elem); }

        inline void setConstant(const value_type &elem) {
            for (size_type i = 0; i < this->size(); ++i) {
                *(begin() + i) = elem;
            }
        }

        inline void setZero() {
            for (size_type i = 0; i < this->size(); ++i) {
                *(begin() + i) = Zero;
            }
        }

        /**
         * \brief Return the number of elements in the Array.
         */
        hur_nodiscard inline size_type size() const noexcept { return mData_.size(); }
        hur_nodiscard inline size_type capacity() const noexcept { return mData_.capacity(); }
        hur_nodiscard inline size_type maxSize() const noexcept { return mData_.maxSize(); }
        hur_nodiscard inline bool empty() const noexcept { return mData_.empty(); }

        hur_nodiscard std::streamsize byteSize() const noexcept {
            return this->size() * sizeof(value_type);
        }

        hur_nodiscard inline const value_type *cdata() const noexcept { return mData_.data(); }
        hur_nodiscard inline const value_type *data() const noexcept { return mData_.data(); }
        hur_nodiscard inline value_type *data() noexcept { return mData_.data(); }

        hur_nodiscard inline value_type &first() noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("first() called on empty ArrayDyn<>");
            }
#endif // HUR_DEBUG
            return mData_[0];
        }

        hur_nodiscard inline const value_type &first() const noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("first() called on empty ArrayDyn<>");
            }
#endif // HUR_DEBUG
            return data()[0];
        }
        hur_nodiscard inline value_type &front() noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("front() called on empty ArrayDyn<>");
            }
#endif // HUR_DEBUG
            return data()[0];
        }
        hur_nodiscard inline const value_type &front() const noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("front() called on empty ArrayDyn<>");
            }
#endif // HUR_DEBUG
            return data()[0];
        }
        hur_nodiscard inline value_type &last() noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("last() called on empty ArrayDyn<>");
            }
#endif // HUR_DEBUG
            return data()[size() - 1];
        }
        hur_nodiscard inline const value_type &last() const noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("last() called on empty ArrayDyn<>");
            }
#endif // HUR_DEBUG
            return data()[size() - 1];
        }
        hur_nodiscard inline value_type &back() noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("back() called on empty ArrayDyn<>");
            }
#endif // HUR_DEBUG
            return data()[size() - 1];
        }
        hur_nodiscard inline const value_type &back() const noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("back() called on empty ArrayDyn<>");
            }
#endif // HUR_DEBUG
            return data()[size() - 1];
        }

        hur_nodiscard inline iterator begin() noexcept { return mData_.data(); }
        hur_nodiscard inline iterator end() noexcept { return (mData_.data() + size()); }

        hur_nodiscard inline const_iterator cbegin() const noexcept { return mData_.data(); }
        hur_nodiscard inline const_iterator cend() const noexcept {
            return (mData_.data() + size());
        }

        hur_nodiscard inline const_iterator begin() const noexcept { return mData_.data(); }
        hur_nodiscard inline const_iterator end() const noexcept {
            return (mData_.data() + size());
        }

        hur_nodiscard inline reverse_iterator rbegin() noexcept {
            return (mData_.data() + size() - 1);
        }
        hur_nodiscard inline reverse_iterator rend() noexcept { return mData_.data() - 1; }

        hur_nodiscard inline const_reverse_iterator crbegin() const noexcept {
            return (mData_.data() + size() - 1);
        }
        hur_nodiscard inline const_reverse_iterator crend() const noexcept {
            return mData_.data() - 1;
        }

        hur_nodiscard inline const_reverse_iterator rbegin() const noexcept {
            return (mData_.data() + size() - 1);
        }
        hur_nodiscard inline const_reverse_iterator rend() const noexcept {
            return mData_.data() - 1;
        }

        hur_nodiscard inline value_type &operator[](const size_type i) noexcept {
            return mData_.element(i);
        }

        hur_nodiscard inline const value_type &operator[](const size_type i) const noexcept {
            return mData_.element(i);
        }

        /**\brief Return element of List.*/
        hur_nodiscard inline value_type &operator()(const size_type i) noexcept {
            return mData_.element(i);
        }

        /**
         *\brief Return element of constant List
         */
        hur_nodiscard inline const value_type &operator()(const size_type i) const noexcept {
            return mData_.element(i);
        }

        hur_nodiscard bool operator==(const ArrayDyn &other) const noexcept {
            if (this->size() != other.size()) {
                return false;
            }
            for (size_type i = 0; i < this->size(); ++i) {
                if (*(begin() + i) != *(other.begin() + i)) {
                    return false;
                }
            }
            return true;
        }

        hur_nodiscard inline bool operator!=(const ArrayDyn &other) const noexcept {
            return !operator==(other);
        }
    }; // namespace OpenHurricane

    template <class Type, typename sizeType> void sort(ArrayDyn<Type, sizeType> &arr) {
        std::sort(arr.begin(), arr.end());
    }
    template <class Type, typename sizeType, class compare>
    void sort(ArrayDyn<Type, sizeType> &arr, const compare &cmp) {
        std::sort(arr.begin(), arr.end(), cmp);
    }

    template <class Type, typename sizeType> void stableSort(ArrayDyn<Type, sizeType> &arr) {
        std::stable_sort(arr.begin(), arr.end());
    }

    template <class Type, typename sizeType, class compare>
    void stableSort(ArrayDyn<Type, sizeType> &ar, const compare &cmp) {
        std::stable_sort(ar.begin(), ar.end(), cmp);
    }

    template <class Type, typename sizeType>
    hur_nodiscard inline std::string toString(const ArrayDyn<Type, sizeType> &arr) {
        std::string str;
        str = "";
        integer i = 0;
        for (const auto &e : arr) {
            str += toString(e);
            if (i < arr.size() - 1) {
                str += ", ";
            }
            i++;
            if ((i % 10 == 0)) {
                str.append("\n");
            }
        }
        return str;
    }
} // namespace OpenHurricane
