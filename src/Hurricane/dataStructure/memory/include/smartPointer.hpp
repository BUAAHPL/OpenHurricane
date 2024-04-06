/*!
 * \file smartPointer.hpp
 * \brief Header of smart pointer.
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
#include "errorAbort.hpp"
#include "integer.hpp"
#include <iostream>

namespace OpenHurricane {
    /**
     * \brief The class of non-copyable pointer to an object.
     */
    template <class Type> class uniquePtr {
    public:
        using pointer = Type *;
        using element_type = Type;

    private:
        pointer ptr_;

    public:
        constexpr inline uniquePtr() noexcept : ptr_(nullptr) {}
        constexpr inline uniquePtr(std::nullptr_t) noexcept : ptr_(nullptr) {}
        inline uniquePtr &operator=(std::nullptr_t) noexcept {
            reset(nullptr);
            return *this;
        }

        inline explicit uniquePtr(pointer ptr) noexcept : ptr_(ptr) {}

        uniquePtr(const uniquePtr &) = delete;
        uniquePtr &operator=(const uniquePtr &) = delete;

        inline uniquePtr(uniquePtr &&other) noexcept : ptr_(other.release()) {}
        inline uniquePtr &operator=(uniquePtr &&other) noexcept {
            reset(other.release());
            return *this;
        }

        inline void swap(uniquePtr &other) noexcept {
            pointer p = ptr_;
            ptr_ = other.ptr_;
            other.ptr_ = p;
        }

        inline ~uniquePtr() noexcept { HurDelete(ptr_); }

        hur_nodiscard inline pointer release() noexcept {
            pointer p = ptr_;
            ptr_ = nullptr;
            return p;
        }

        inline pointer unsafeRelease() noexcept {
            pointer p = ptr_;
            ptr_ = nullptr;
            return p;
        }

        hur_nodiscard inline pointer setRelease(std::nullptr_t) noexcept { return release(); }
        hur_nodiscard inline pointer setRelease(pointer ptr) noexcept {
            pointer old = ptr_;
            ptr_ = ptr;
            return old;
        }
        inline void reset(std::nullptr_t) noexcept { HurDelete(ptr_); }
        inline void reset(pointer ptr) noexcept {
            pointer old = ptr_;
            ptr_ = ptr;
            HurDelete(old);
        }

        inline void clear() noexcept { HurDelete(ptr_); }

        hur_nodiscard inline pointer operator->() const noexcept { return ptr_; }
        hur_nodiscard inline const pointer get() const noexcept { return ptr_; }
        hur_nodiscard inline pointer get() noexcept { return ptr_; }
        hur_nodiscard inline const element_type &operator*() const noexcept { return *ptr_; }

        hur_nodiscard inline element_type &operator*() noexcept { return *ptr_; }

        inline explicit operator bool() const noexcept { return ptr_ != nullptr; }

        hur_nodiscard inline bool isNull() const noexcept { return ptr_ == nullptr; }
    };

    /*!\brief Class for reference counted pointer management.*/
    template <class Type> class sharedPtr {
    public:
        using pointer = Type *;
        using element_type = Type;
        using count_type = unsigned int;

    private:
        pointer ptr_;
        count_type *sharedCount_;

    public:
        constexpr inline sharedPtr() noexcept : ptr_(nullptr), sharedCount_(nullptr) {}

        constexpr inline sharedPtr(std::nullptr_t) noexcept
            : ptr_(nullptr), sharedCount_(nullptr) {}

        inline explicit sharedPtr(pointer ptr) noexcept
            : ptr_(ptr), sharedCount_(new count_type(1)) {}

        inline sharedPtr(const sharedPtr &other) noexcept {
            ptr_ = other.ptr_;
            sharedCount_ = other.sharedCount_;
            if (sharedCount_ != nullptr) {
                (*sharedCount_)++;
            }
        }

        inline sharedPtr &operator=(const sharedPtr &other) noexcept {
            if (ptr_ != other.ptr_) {
                if (sharedCount_ != nullptr) {
                    (*sharedCount_)--;
                    if (*sharedCount_ == 0) {
                        HurDelete(ptr_);
                        HurDelete(sharedCount_);
                    }
                }

                ptr_ = other.ptr_;
                sharedCount_ = other.sharedCount_;
                if (sharedCount_) {
                    (*sharedCount_)++;
                }
            }
            return *this;
        }

        inline sharedPtr(sharedPtr &&other) noexcept : ptr_(nullptr), sharedCount_(nullptr) {
            this->swap(other);
        }

        inline sharedPtr &operator=(sharedPtr &&other) noexcept {
            sharedPtr(std::move(other)).swap(*this);
        }

        inline ~sharedPtr() noexcept {
            if (sharedCount_ != nullptr) {
                (*sharedCount_)--;
                if (*sharedCount_ == 0) {
                    HurDelete(ptr_);
                    HurDelete(sharedCount_);
                }
            }
        }

        inline void swap(sharedPtr &other) noexcept {
            std::swap(ptr_, other.ptr_);
            std::swap(sharedCount_, other.sharedCount_);
        }

        hur_nodiscard inline pointer operator->() const noexcept { return ptr_; }
        hur_nodiscard inline const pointer get() const noexcept { return ptr_; }
        hur_nodiscard inline pointer get() noexcept { return ptr_; }
        hur_nodiscard inline const element_type &operator*() const noexcept { return *ptr_; }
        hur_nodiscard inline element_type &operator*() noexcept { return *ptr_; }

        inline explicit operator bool() const noexcept { return ptr_ != nullptr; }

        hur_nodiscard inline bool isNull() const noexcept { return ptr_ == nullptr; }

        hur_nodiscard inline count_type sharedCount() const noexcept {
            return sharedCount_ != nullptr ? *sharedCount_ : 0;
        }

        inline void reset(std::nullptr_t) noexcept { sharedPtr().swap(*this); }
        inline void reset(pointer ptr) noexcept { sharedPtr(ptr).swap(*this); }

        inline pointer unsafeRelease() noexcept {
            if (sharedCount_ != nullptr) {
                (*sharedCount_)--;
                if (*sharedCount_ == 0) {
                    HurDelete(sharedCount_);
                } else {
                    sharedCount_ = nullptr;
                }
            }
            auto tmp = ptr_;
            ptr_ = nullptr;
            return tmp;
        }
    };

} // namespace OpenHurricane