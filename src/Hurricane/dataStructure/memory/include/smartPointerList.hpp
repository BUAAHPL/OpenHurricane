/*!
 * \file smartPointerList.hpp
 * \brief Header of list of smart pointer.
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
#include "List.hpp"
#include "smartPointer.hpp"

namespace OpenHurricane {
    template <class Type> class plainPointerList;
    template <class Type> class uniquePtrList;
    template <class Type> class sharedPtrList;

    template <class Type> using PtrList = plainPointerList<Type>;
    template <class Type> using uPtrList = uniquePtrList<Type>;
    template <class Type> using sPtrList = sharedPtrList<Type>;

    template <class Type> class plainPointerList {
    public:
        using Base = List<Type *>;
        using value_type = Type;
        using reference = Type &;
        using const_reference = const Type &;
        using difference_type = typename Base::difference_type;
        using size_type = typename Base::size_type;
        using pointer = Type *;
        using element_type = Type;

    private:
        Base ptr_;

    public:
        inline plainPointerList() : ptr_() {}
        inline explicit plainPointerList(const size_type s) : ptr_(s) {
            for (size_type i = 0; i < s; ++i) {
                ptr_[i] = nullptr;
            }
        }
        inline plainPointerList(const plainPointerList &other) : ptr_(other.size()) {
            for (size_type i = 0; i < other.size(); ++i) {
                if (other.ptr_[i] == nullptr) {
                    ptr_[i] = nullptr;
                } else {
                    ptr_[i] = other.ptr_[i]->clone().release();
                }
            }
        }
        plainPointerList &operator=(const plainPointerList &other) {
            if (this != std::addressof(other)) {
                if (!empty()) {
                    resize(other.size());
                    for (integer i = 0; i < size(); ++i) {
                        ptr_[i] = other[i].clone().release();
                    }
                } else if (other.size() == this->size()) {
                    for (integer i = 0; i < size(); ++i) {
                        if (ptr_[i] != nullptr) {
                            (*this)[i] = other[i];
                        } else {
                            ptr_[i] = other[i].clone().release();
                        }
                    }
                } else {
                    clear();
                    resize(other.size());
                    for (integer i = 0; i < size(); ++i) {
                        ptr_[i] = other[i].clone().release();
                    }
                }
            }
            return *this;
        }
        template <class CloneArgs>
        inline plainPointerList(const plainPointerList &other, const CloneArgs &cloneArgs)
            : ptr_(other.size()) {
            for (size_type i = 0; i < other.size(); ++i) {
                if (other.ptr_[i] == nullptr) {
                    ptr_[i] = nullptr;
                } else {
                    ptr_[i] = other.ptr_[i]->clone(cloneArgs).release();
                }
            }
        }
        inline plainPointerList(plainPointerList &&other) noexcept : ptr_(std::move(other.ptr_)) {}
        plainPointerList &operator=(plainPointerList &&other) noexcept { transfer(other); }
        inline ~plainPointerList() noexcept { clear(); }

        hur_nodiscard inline size_type size() const noexcept { return ptr_.size(); }
        hur_nodiscard inline size_type capacity() const noexcept { return ptr_.capacity(); }
        hur_nodiscard inline size_type maxSize() const noexcept { return ptr_.maxSize(); }
        hur_nodiscard inline bool empty() const noexcept { return ptr_.empty(); }
        hur_nodiscard inline reference first() noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("first() called on empty plainPointerList<>");
            }
#endif // HUR_DEBUG
            return *(ptr_.first());
        }

        hur_nodiscard inline const_reference first() const noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("first() called on empty plainPointerList<>");
            }
#endif // HUR_DEBUG
            return *(ptr_.first());
        }

        /**\brief Return the first element of the list.*/
        hur_nodiscard inline reference front() noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("front() called on empty plainPointerList<>");
            }
#endif // HUR_DEBUG
            return *(ptr_.first());
        }

        /**\brief Return first element of the list.*/
        hur_nodiscard inline const_reference front() const noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("front() called on empty plainPointerList<>");
            }
#endif // HUR_DEBUG
            return *(ptr_.first());
        }

        /**\brief Return the last element of the list.*/
        hur_nodiscard inline reference last() noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("last() called on empty plainPointerList<>");
            }
#endif // HUR_DEBUG
            return *(ptr_.last());
        }

        /**\brief Return the last element of the list.*/
        hur_nodiscard inline const_reference last() const noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("last() called on empty plainPointerList<>");
            }
#endif // HUR_DEBUG
            return *(ptr_.last());
        }

        /**\brief Return the last element of the list.*/
        hur_nodiscard inline reference back() noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("last() called on empty plainPointerList<>");
            }
#endif // HUR_DEBUG
            return *(ptr_.last());
        }

        /**\brief Return the last element of the list.*/
        hur_nodiscard inline const_reference back() const noexcept {
#ifdef HUR_DEBUG
            if (empty()) {
                LFatal("last() called on empty plainPointerList<>");
            }
#endif // HUR_DEBUG
            return *(ptr_.last());
        }

        void resize(const size_type newSize) {
            integer oldSize = size();
            if (newSize <= 0) {
                clear();
            } else if (newSize < oldSize) {
                for (integer i = newSize; i < oldSize; ++i) {
                    HurDelete(ptr_[i]);
                }
                ptr_.resize(newSize);
            } else if (newSize > oldSize) {
                ptr_.resize(newSize);
                for (integer i = oldSize; i < newSize; ++i) {
                    ptr_[i] = nullptr;
                }
            }
        }
        void clear() noexcept {
            for (integer i = 0; i < this->size(); i++) {
                HurDelete(ptr_[i]);
            }
            ptr_.clear();
        }
        void unbind() noexcept {
            for (integer i = 0; i < this->size(); i++) {
                if (ptr_[i] != nullptr) {
                    ptr_[i] = nullptr;
                }
            }
        }

        void swap(plainPointerList &other) noexcept { ptr_.swap(other.ptr_); }
        void transfer(plainPointerList &other) noexcept {
            clear();
            swap(other);
        }
        hur_nodiscard inline bool set(const size_type i) const { return ptr_[i] != nullptr; }
        inline pointer set(const size_type i, pointer newElePtr) {
            pointer oldPtr = ptr_[i];
            ptr_[i] = newElePtr;
            return oldPtr;
        }
        inline void append(pointer newElePtr) {
            integer sz = size();
            resize(sz + 1);
            ptr_[sz] = newElePtr;
        }

        inline void append(uniquePtr<Type> &&aptr) { append(aptr.release()); }
        hur_nodiscard inline reference operator[](const size_type i) {
#ifdef HUR_DEBUG
            if (ptr_[i] == nullptr) {
                LFatal("Hanging pointer(null pointer) at index: %d", i);
            }
#endif // HUR_DEBUG
            return *(ptr_[i]);
        }
        hur_nodiscard inline const_reference operator[](const size_type i) const {
#ifdef HUR_DEBUG
            if (ptr_[i] == nullptr) {
                LFatal("Hanging pointer(null pointer) at index: %d", i);
            }
#endif // HUR_DEBUG
            return *(ptr_[i]);
        }

        hur_nodiscard inline pointer operator()(const size_type i) { return ptr_[i]; }
        hur_nodiscard inline const pointer operator()(const size_type i) const { return ptr_[i]; }

        class iterator;
        friend class iterator;

        class iterator {
        private:
            Type **pointer_;

        public:
            inline iterator(Type **ptr) : pointer_(ptr) {}

            hur_nodiscard inline bool operator==(const iterator &iter) const noexcept {
                return pointer_ == iter.pointer_;
            }
            hur_nodiscard inline bool operator!=(const iterator &iter) const noexcept {
                return pointer_ != iter.pointer_;
            }

            hur_nodiscard inline Type &operator*() noexcept { return **pointer_; }
            hur_nodiscard inline Type &operator()() noexcept { return operator*(); }

            inline iterator operator++() noexcept {
                ++pointer_;
                return *this;
            }
            inline iterator operator++(int) noexcept {
                iterator tmp = *this;
                ++pointer_;
                return tmp;
            }

            inline iterator operator--() noexcept {
                --pointer_;
                return *this;
            }

            inline iterator operator--(int) noexcept {
                iterator tmp = *this;
                --pointer_;
                return tmp;
            }

            inline iterator operator+=(size_type n) noexcept {
                pointer_ += n;
                return *this;
            }

            inline iterator operator-=(size_type n) noexcept {
                pointer_ -= n;
                return *this;
            }

            hur_nodiscard inline Type &operator[](size_type n) noexcept { return *(*this + n); }

            hur_nodiscard inline bool operator<(const iterator &iter) const noexcept {
                return pointer_ < iter.pointer_;
            }
            hur_nodiscard inline bool operator>(const iterator &iter) const noexcept {
                return pointer_ > iter.pointer_;
            }

            hur_nodiscard inline bool operator<=(const iterator &iter) const noexcept {
                return pointer_ <= iter.pointer_;
            }
            hur_nodiscard inline bool operator>=(const iterator &iter) const noexcept {
                return pointer_ >= iter.pointer_;
            }
        };

        hur_nodiscard inline iterator begin() noexcept { return ptr_.begin(); }
        hur_nodiscard inline iterator end() noexcept { return ptr_.end(); }

        class const_iterator;
        friend class const_iterator;

        class const_iterator {
        private:
            const Type *const *pointer_;

        public:
            inline const_iterator(const Type *const *ptr) : pointer_(ptr) {}

            inline const_iterator(const iterator &iter) : pointer_(iter.pointer_) {}

            hur_nodiscard inline bool operator==(const const_iterator &iter) const noexcept {
                return pointer_ == iter.pointer_;
            }

            hur_nodiscard inline bool operator!=(const const_iterator &iter) const noexcept {
                return pointer_ != iter.pointer_;
            }

            hur_nodiscard inline const Type &operator*() noexcept { return **pointer_; }
            hur_nodiscard inline const Type &operator()() noexcept { return operator*(); }

            inline const_iterator operator++() noexcept {
                ++pointer_;
                return *this;
            }
            inline const_iterator operator++(int) noexcept {
                const_iterator tmp = *this;
                ++pointer_;
                return tmp;
            }

            inline const_iterator operator--() noexcept {
                --pointer_;
                return *this;
            }
            inline const_iterator operator--(int) noexcept {
                const_iterator tmp = *this;
                --pointer_;
                return tmp;
            }

            inline const_iterator operator+=(size_type n) noexcept {
                pointer_ += n;
                return *this;
            }

            inline const_iterator operator-=(size_type n) noexcept {
                pointer_ -= n;
                return *this;
            }

            hur_nodiscard inline Type &operator[](size_type n) noexcept { return *(*this + n); }

            hur_nodiscard inline bool operator<(const const_iterator &iter) const noexcept {
                return pointer_ < iter.pointer_;
            }
            hur_nodiscard inline bool operator>(const const_iterator &iter) const noexcept {
                return pointer_ > iter.pointer_;
            }

            hur_nodiscard inline bool operator<=(const const_iterator &iter) const noexcept {
                return pointer_ <= iter.pointer_;
            }
            hur_nodiscard inline bool operator>=(const const_iterator &iter) const noexcept {
                return pointer_ >= iter.pointer_;
            }
        };

        hur_nodiscard inline const_iterator begin() const noexcept { return ptr_.begin(); }

        hur_nodiscard inline const_iterator end() const noexcept { return ptr_.end(); }

        hur_nodiscard inline const_iterator cbegin() const noexcept { return ptr_.cbegin(); }

        hur_nodiscard inline const_iterator cend() const noexcept { return ptr_.cend(); }
    };

    template <class Type> class uniquePtrList {
    public:
        using Base = List<uniquePtr<Type>>;
        using size_type = typename Base::size_type;
        using value_type = typename Base::value_type;

        using pointer = Type *;
        using element_type = Type;

    private:
        uniquePtr<Base> ptr_;

    public:
        inline uniquePtrList() : ptr_() {}

        inline explicit uniquePtrList(const size_type s) : ptr_(new Base(s)) {}

        uniquePtrList(const uniquePtrList &) = delete;
        uniquePtrList &operator=(const uniquePtrList &) = delete;

        inline uniquePtrList(uniquePtrList &&other) noexcept : ptr_(std::move(other.ptr_)) {}
        inline uniquePtrList &operator=(uniquePtrList &&other) noexcept {
            ptr_.clear();
            ptr_.swap(other.ptr_);
            return *this;
        }

        inline ~uniquePtrList() noexcept {}

        hur_nodiscard inline size_type size() const noexcept {
            if (ptr_) {
                return ptr_->size();
            } else {
                return size_type(0);
            }
        }
        hur_nodiscard inline bool empty() const noexcept {
            if (ptr_) {
                return ptr_->empty();
            } else {
                return true;
            }
        }

        inline void clear() noexcept {
            if (ptr_) {
                ptr_->clear();
            }
        }
        inline void unbind() {
            if (ptr_) {
                for (integer i = 0; i < ptr_->size(); ++i) {
                    pointer tmpPtr = (*ptr_)[i].release();
                }
            }
        }
        inline void swap(uniquePtrList &other) noexcept { other.ptr_.swap(this->ptr_); }
        inline void transfer(uniquePtrList &other) noexcept {
            clear();
            other.ptr_.transfer(this->ptr_);
        }
        void resize(const size_type newSize) {
            clear();
            if (newSize < 0) {
                LFatal("invalid index");
            }
            auto oldSize = size();
            if (newSize > oldSize) {
                ptr_.clear();
                ptr_ = uniquePtr<Base>(new Base(newSize));
            }
        }

        hur_nodiscard inline bool set(const size_type i) const {
            if (ptr_) {
                return (*ptr_)[i];
            } else {
                return false;
            }
        }
        inline value_type set(const size_type i, pointer pi) {
            if (!ptr_) {
                ptr_ = uniquePtr<Base>(new Base(i + 1));
            }
            value_type tmpP(pi);
            tmpP.swap((*ptr_)[i]);
            return tmpP;
        }

        inline value_type set(const size_type i, value_type &&pi) {
            if (!ptr_) {
                ptr_ = uniquePtr<Base>(new Base(i + 1));
            }
            value_type tmpP(std::move(pi));
            tmpP.swap((*ptr_)[i]);
            return tmpP;
        }

        hur_nodiscard inline element_type &operator[](const size_type i) { return *(*ptr_)[i]; }

        hur_nodiscard inline const element_type &operator[](const size_type i) const {
            return *(*ptr_)[i];
        }

        hur_nodiscard inline pointer operator()(const size_type i) { return (*ptr_)[i].get(); }

        hur_nodiscard inline const pointer operator()(const size_type i) const {
            { return (*ptr_)[i].get(); }
        }
    };

    template <class Type> class sharedPtrList {
    public:
        using Base = List<sharedPtr<Type>>;
        using size_type = typename Base::size_type;
        using value_type = typename Base::value_type;

        using pointer = Type *;
        using element_type = Type;

    private:
        Base ptr_;

    public:
        inline sharedPtrList() : ptr_() {}

        inline explicit sharedPtrList(const size_type s) : ptr_(s) {}
        inline sharedPtrList(const sharedPtrList &other) : ptr_(other.ptr_) {}
        inline sharedPtrList &operator=(const sharedPtrList &other) {
            if (this != std::addressof(other)) {
                ptr_ = other.ptr_;
            }
            return *this;
        }

        inline sharedPtrList(sharedPtrList &&other) noexcept : ptr_(std::move(other.ptr_)) {}
        inline sharedPtrList &operator=(sharedPtrList &&other) noexcept {
            ptr_.clear();
            ptr_.swap(other.ptr_);
            return *this;
        }
        inline ~sharedPtrList() noexcept {}

        hur_nodiscard inline size_type size() const noexcept { return ptr_.size(); }
        hur_nodiscard inline bool empty() const noexcept { return ptr_.empty(); }

        inline void clear() noexcept { ptr_.clear(); }
        inline void swap(sharedPtrList &other) noexcept { other.ptr_.swap(this->ptr_); }
        inline void transfer(sharedPtrList &other) noexcept {
            clear();
            other.ptr_.transfer(this->ptr_);
        }

        void resize(const size_type newSize) {
            if (newSize < 0) {
                LFatal("invalid index");
            }
            auto oldSize = size();
            if (newSize <= 0) {
                clear();
            } else if (newSize < oldSize) {
                for (integer i = newSize; i < oldSize; ++i) {
                    ptr_[i].reset(nullptr);
                }
                ptr_.resize(newSize);
            } else if (newSize > oldSize) {
                ptr_.resize(newSize);
            }
        }

        hur_nodiscard inline bool set(const size_type i) const { return ptr_[i]; }
        inline void set(const size_type i, pointer pi) { ptr_[i].reset(pi); }
        inline void set(const size_type i, const value_type &pi) { ptr_[i] = pi; }
        inline void set(const size_type i, value_type &&pi) { ptr_[i].swap(pi); }
        inline void append(pointer pi) { ptr_.append(sharedPtr<Type>(pi)); }
        inline void append(const value_type &pi) { ptr_.append(pi); }
        inline void append(uniquePtr<Type> &&pi) { ptr_.append(sharedPtr<Type>(pi.release())); }

        hur_nodiscard inline element_type &operator[](const size_type i) { return *ptr_[i]; }

        hur_nodiscard inline const element_type &operator[](const size_type i) const {
            return *ptr_[i];
        }

        hur_nodiscard inline pointer operator()(const size_type i) { return ptr_[i].get(); }

        hur_nodiscard inline const pointer operator()(const size_type i) const {
            { return ptr_[i].get(); }
        }

        inline void unbind() {
            if (!ptr_.empty()) {
                for (integer i = 0; i < ptr_.size(); ++i) {
                    ptr_[i].unsafeRelease();
                }
            }
        }
    };
} // namespace OpenHurricane
