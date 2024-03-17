/*!
 * \file ArrayStorage.hpp
 * \brief Headers of the class of Array storage.
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
#include "errorAbort.hpp"
#include "integer.hpp"
#include "zero.hpp"
#include <iostream>

namespace OpenHurricane {

    /**
     * \brief Purely fix-size array.
     */
    template <class Type, int sizeOfElement> class storageFix {
    public:
        using value_type = Type;
        using reference = Type &;
        using const_reference = const Type &;
        using difference_type = int;
        using size_type = int;

        static constexpr int nElements_ = sizeOfElement;

    protected:
        static_assert(sizeOfElement > 0 && sizeOfElement <= INT_MAX,
                      "Size must be positive (>0) and also fit as a signed value");

        Type element_[sizeOfElement];

    public:
        /**
         * \brief Null constructor.
         */
        constexpr inline storageFix() : element_() {}

        constexpr inline storageFix(const value_type &elem) {
            for (int i = 0; i < sizeOfElement; ++i) {
                element_[i] = elem;
            }
        }

        constexpr inline storageFix(const zero) : element_() {
            for (int i = 0; i < sizeOfElement; ++i) {
                element_[i] = Zero;
            }
        }

        /**
         * \brief Copy constructor.
         */
        constexpr inline storageFix(const storageFix &other) {
            for (int i = 0; i < sizeOfElement; ++i) {
                element_[i] = other.element_[i];
            }
        }

        constexpr inline storageFix &operator=(const storageFix &other) {
            if (this != std::addressof(other)) {
                for (int i = 0; i < sizeOfElement; ++i) {
                    element_[i] = other.element_[i];
                }
            }
            return *this;
        }

        template <class otherType>
        inline explicit storageFix(const storageFix<otherType, sizeOfElement> &other) {
            for (int i = 0; i < sizeOfElement; ++i) {
                element_[i] = static_cast<const Type &>(other.element_[i]);
            }
        }

        inline ~storageFix() noexcept {}

        void swap(storageFix &other) noexcept { std::swap(element_, other.element_); }

        void clear() noexcept {}

        /**
         * \brief Remove all the entries that are same with value of the Array.
         */
        constexpr inline void remove(const value_type &elem) {
            for (size_type i = 0; i < sizeOfElement; i++) {
                if (element_[i] == elem) {
                    element_[i] = Zero;
                }
            }
        }
        /**
         * \brief Return the number of elements.
         */
        hur_nodiscard inline static constexpr size_type size() noexcept { return sizeOfElement; }

        hur_nodiscard constexpr inline const value_type *data() const noexcept {
            return &element_[0];
        }

        hur_nodiscard constexpr inline value_type *data() noexcept { return &element_[0]; }

        hur_nodiscard constexpr inline const value_type &element(const size_type i) const noexcept {
#ifdef HUR_DEBUG
            if (i >= sizeOfElement || i < 0) {
                LFatal("Index out of range 0~%d while i = %d", sizeOfElement - 1, i);
            }
#endif // HUR_DEBUG
            return element_[i];
        }

        hur_nodiscard constexpr inline value_type &element(const size_type i) noexcept {
#ifdef HUR_DEBUG
            if (i >= sizeOfElement || i < 0) {
                LFatal("Index out of range 0~%d while i = %d", sizeOfElement - 1, i);
            }
#endif // HUR_DEBUG
            return element_[i];
        }

        constexpr inline void element(value_type &ele, const size_type i) const noexcept {
#ifdef HUR_DEBUG
            if (i >= sizeOfElement || i < 0) {
                LFatal("Index out of range 0~%d while i = %d", sizeOfElement - 1, i);
            }
#endif // HUR_DEBUG
            ele = element_[i];
        }

        constexpr inline void replace(const size_type i, const value_type &ele) noexcept {
#ifdef HUR_DEBUG
            if (i >= sizeOfElement || i < 0) {
                LFatal("Index out of range 0~%d while i = %d", sizeOfElement - 1, i);
            }
#endif // HUR_DEBUG
            element_[i] = ele;
        }

        hur_nodiscard constexpr inline const value_type &
        operator[](const size_type i) const noexcept {
            return element(i);
        }
        hur_nodiscard constexpr inline value_type &operator[](const size_type i) noexcept {
            return element(i);
        }
    };

    /**
     * \brief Purely dynamic array.
     */
    template <class Type, typename sizeType> class storageDyn {
    public:
        using value_type = Type;
        using reference = Type &;
        using const_reference = const Type &;
        using difference_type = sizeType;
        using size_type = sizeType;

    private:
        size_type size_;

        value_type *hur_restrict element_;

        size_type capacity_;

        inline void allocate() {
            HurDeleteDynArray(element_);
            if (capacity_ > 0) {
                element_ = new value_type[capacity_];
            } else {
                capacity_ = 0;
                size_ = 0;
            }
        }

        inline void checkSize(const size_type newSize) const {
            if (newSize < 0) {
                errorAbortStr(("Bad size: " + toString(newSize)));
            }
        }

    public:
        /**
         * \brief Null constructor.
         */
        inline storageDyn() : size_(0), element_(nullptr), capacity_(0) {}

        inline storageDyn(const size_type size) : size_(size), element_(nullptr), capacity_(size) {
            checkSize(size);
            allocate();
        }

        inline storageDyn(const size_type size, const size_type capacity)
            : size_(size), element_(nullptr), capacity_(capacity) {
            checkSize(size);
            checkSize(capacity);
            allocate();
        }
        inline storageDyn(const size_type &s, const size_type &s2, const value_type &elem)
            : size_(s), element_(nullptr), capacity_(s2) {
            checkSize(s);
            checkSize(s2);
            allocate();
            for (int i = 0; i < size_; ++i) {
                element_[i] = elem;
            }
        }

        inline storageDyn(const size_type &s, const size_type &s2, const zero)
            : size_(s), element_(nullptr), capacity_(s2) {
            checkSize(s);
            checkSize(s2);
            allocate();
            for (int i = 0; i < size_; ++i) {
                element_[i] = Zero;
            }
        }

        /**
         * \brief Copy constructor.
         */
        inline storageDyn(const storageDyn &other)
            : size_(other.size_), element_(nullptr), capacity_(other.capacity_) {
            allocate();
            if (element_ != nullptr) {
                for (size_type i = 0; i < size_; ++i) {
                    *(element_ + i) = *(other.element_ + i);
                }
            }
        }

        inline storageDyn &operator=(const storageDyn &other) {
            if (this != std::addressof(other)) {
                size_ = other.size_;
                if (size_ > capacity_) {
                    capacity_ = other.capacity_;
                    allocate();
                }
                for (size_type i = 0; i < size_; ++i) {
                    *(element_ + i) = *(other.element_ + i);
                }
            }
            return *this;
        }

        template <class otherType>
        inline explicit storageDyn(const storageDyn<otherType, sizeType> &other)
            : size_(other.size_), element_(nullptr), capacity_(other.capacity_) {
            for (size_type i = 0; i < size_; ++i) {
                *(element_ + i) = static_cast<otherType>(*(other.element_ + i));
            }
        }

        template <class otherType>
        inline storageDyn &operator=(const storageDyn<otherType, sizeType> &other) {
            if (this != std::addressof(other)) {
                size_ = other.size_;
                if (size_ > capacity_) {
                    capacity_ = other.capacity_;
                    allocate();
                }
                for (size_type i = 0; i < size_; ++i) {
                    *(element_ + i) = static_cast<otherType>(*(other.element_ + i));
                }
            }
            return *this;
        }

        /**
         * \brief Construct as copy from rvalue.
         */
        inline storageDyn(storageDyn &&other) noexcept
            : size_(std::move(other.size_)), element_(std::move(other.element_)),
              capacity_(std::move(other.capacity_)) {
            other.size_ = 0;
            other.element_ = nullptr;
            other.capacity_ = 0;
        }
        inline storageDyn &operator=(storageDyn &&other) noexcept {
            std::swap(size_, other.size_);
            std::swap(element_, other.element_);
            std::swap(capacity_, other.capacity_);
            return *this;
        }

        inline void swap(storageDyn &other) noexcept {
            std::swap(size_, other.size_);
            std::swap(element_, other.element_);
            std::swap(capacity_, other.capacity_);
        }

        inline void transfer(storageDyn &other) noexcept {
            this->swap(other);
            other.clear();
        }

        inline ~storageDyn() noexcept { clear(); }

        inline void clear() noexcept {
            size_ = 0;
            capacity_ = 0;
            HurDeleteDynArray(element_);
        }

        hur_nodiscard inline size_type size() const noexcept { return size_; }
        hur_nodiscard inline size_type capacity() const noexcept { return capacity_; }
        hur_nodiscard inline size_type maxSize() const noexcept { return integerMax; }
        hur_nodiscard inline bool empty() const noexcept { return size_ <= 0; }

    private:
        inline void checkIndex(const size_type i) const noexcept {
#ifdef HUR_DEBUG
            if (size_ <= 0) {
                LFatal("Attempt to access element from empty array");
            }
            if (i >= size_ || i < 0) {
                LFatal("Index out of range 0~%d while i = %d", size_ - 1, i);
            }
#endif // HUR_DEBUG
        }

    public:
        hur_nodiscard inline const value_type *data() const noexcept { return element_; }
        hur_nodiscard inline value_type *data() noexcept { return element_; }
        hur_nodiscard inline const value_type &element(const size_type i) const noexcept {
            checkIndex(i);
            return element_[i];
        }

        hur_nodiscard inline value_type &element(const size_type i) noexcept {
            checkIndex(i);
            return element_[i];
        }

        inline void element(value_type &ele, const size_type i) const noexcept {
            checkIndex(i);
            ele = element_[i];
        }

        inline void replace(const size_type i, const value_type &ele) noexcept {
            checkIndex(i);
            element_[i] = ele;
        }
        hur_nodiscard inline value_type &operator[](const size_type i) noexcept {
            checkIndex(i);
            return element_[i];
        }

        hur_nodiscard inline const value_type &operator[](const size_type i) const noexcept {
            checkIndex(i);
            return element_[i];
        }

        inline void resize(const size_type newSize) {
            checkSize(newSize);
            if (newSize > capacity_) {
                storageDyn tmpArray(newSize);
                auto minSize = min(newSize, size_);
                for (size_type i = 0; i < minSize; ++i) {
                    *(tmpArray.element_ + i) = *(element_ + i);
                }
                this->swap(tmpArray);
            } else {
                size_ = newSize;
            }
        }

        inline void resize(const size_type newSize, const value_type &elem) {
            checkSize(newSize);
            if (newSize > capacity_) {
                storageDyn tmpArray(newSize);
                this->swap(tmpArray);
            } else {
                size_ = newSize;
            }
            for (size_type i = 0; i < size_; ++i) {
                *(element_ + i) = elem;
            }
        }

        inline void resize(const size_type newSize, const zero) {
            checkSize(newSize);
            if (newSize > capacity_) {
                storageDyn tmpArray(newSize);
                this->swap(tmpArray);
            } else {
                size_ = newSize;
            }
            for (size_type i = 0; i < size_; ++i) {
                *(element_ + i) = Zero;
            }
        }

        inline void reserve(const size_type newCapacity) {
            if (newCapacity > capacity_) {
                storageDyn tmpArray(newCapacity);
                tmpArray.size_ = size_;
                for (size_type i = 0; i < size_; ++i) {
                    *(tmpArray.element_ + i) = *(element_ + i);
                }
                this->swap(tmpArray);
            }
        }
        inline void shrinkToFit() {
            if (size_ != capacity_) {
                storageDyn tmpArray(size_);
                for (size_type i = 0; i < size_; ++i) {
                    *(tmpArray.element_ + i) = *(element_ + i);
                }
                this->swap(tmpArray);
            }
        }
        /**
         * \brief Remove all the entries that are same with value of the Array.
         */
        inline void remove(const value_type &elem) {
            size_type *findIdList = new size_type[size_];
            size_type findNum = 0;
            for (size_type i = 0; i < size_; i++) {
                if (element_[i] == elem) {
                    findIdList[findNum] = i;
                    findNum++;
                }
            }
            if (findNum != 0) {
                size_type newSize = size_ - findNum;
                storageDyn tmpArray(newSize, capacity_);
                size_type j = 0;
                value_type *oldEle = &element_[0];
                value_type *newEle = &tmpArray.element_[0];
                for (size_type i = 0; i < size_; i++) {
                    if (findIdList[j] == i) {
                        oldEle++;
                        j++;
                        if (j > findNum) {
                            break;
                        }
                        continue;
                    }
                    *newEle++ = *oldEle++;
                }
                this->swap(tmpArray);
            }
            HurDeleteDynArray(findIdList);
        }

        /**
         * \brief Insert an element at the entry specified by iter of the array.
         * \param[in] iter - The entry where to insert.
         * \param[in] elem - The value to be inserted.
         */
        void insert(const value_type *iter, const value_type &elem) {
            if ((iter > &element_[size_]) || (iter < &element_[0])) {
                LFatal("Beyond the size, please check!");
            }
            size_type newSize = size_ + 1;
            storageDyn tmpArray(newSize, max(capacity_, newSize));

            value_type *oldElem = &element_[0];
            value_type *newElem = &tmpArray.element_[0];
            size_type i = size_;
            while (i--) {
                if (oldElem == iter) {
                    *newElem++ = elem;
                }
                *newElem++ = *oldElem++;
            }
            if (iter == &element_[size_]) {
                *newElem++ = elem;
            }
            this->swap(tmpArray);
        }

        /**
         * \brief Insert an element at the entry specified by iter of the array.
         * \param[in] iter - The entry where to insert.
         * \param[in] num - Specify the total number to be insert.
         * \param[in] elem - The value to be inserted.
         */
        void insert(const value_type *iter, const size_type num, const value_type &elem) {
            if (num >= 1) {
                if ((iter > &element_[size_]) || (iter < &element_[0])) {
                    LFatal("Beyond the size, please check!");
                }
                size_type newSize = this->size_ + num;
                storageDyn tmpArray(newSize, max(capacity_, newSize));

                value_type *oldElem = &element_[0];
                value_type *newElem = &tmpArray.element_[0];
                size_type i = size_;
                while (i--) {
                    if (oldElem == iter) {
                        for (size_type j = 0; j < num; j++) {
                            *newElem++ = elem;
                        }
                    }
                    *newElem++ = *oldElem++;
                }
                if (iter == &element_[size_]) {
                    for (size_type j = 0; j < num; j++) {
                        *newElem++ = elem;
                    }
                }
                this->swap(tmpArray);
            }
        }
    };

} // namespace OpenHurricane