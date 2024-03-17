/*!
 * \file bounded.hpp
 * \brief Headers of bounded data.
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
#include "smartPointer.hpp"
namespace OpenHurricane {
    /**
     * \brief Class of bounded data.
     */
    template <class Type> class bounded {
    public:
        enum boundedTypes : short { BOUNDED_WITH_VALUE, BOUNDED_FLAG };

    private:
        boundedTypes boundedType_;

        /** \brief Is data lower bounded? */
        bool isLowerBounded_;

        /** \brief Is data upper bounded? */
        bool isUpperBounded_;

        /** \brief The lower bounded value. */
        uniquePtr<Type> lowerBoundedPtr_;

        /** \brief The upper bounded value. */
        uniquePtr<Type> upperBoundedPtr_;

    public:
        /**
         * \brief Null constructor.
         */
        inline bounded()
            : boundedType_(BOUNDED_WITH_VALUE), isLowerBounded_(false), isUpperBounded_(false),
              lowerBoundedPtr_(nullptr), upperBoundedPtr_(nullptr) {}

        /**
         * \brief Copy constructor.
         */
        inline bounded(const bounded<Type> &b)
            : boundedType_(b.boundedType_), isLowerBounded_(b.isLowerBounded_),
              isUpperBounded_(b.isUpperBounded_), lowerBoundedPtr_(nullptr),
              upperBoundedPtr_(nullptr) {
            if (isLowerBounded_) {
                if (!b.lowerBoundedPtr_) {
                    LFatal("The pointer for lower bounded value not set");
                }
                lowerBoundedPtr_.reset(new Type());
                *lowerBoundedPtr_ = *b.lowerBoundedPtr_;
            }
            if (isUpperBounded_) {
                if (!b.upperBoundedPtr_) {
                    LFatal("The pointer for upper bounded value not set");
                }
                upperBoundedPtr_.reset(new Type());
                *upperBoundedPtr_ = *b.upperBoundedPtr_;
            }

            checkBoundValue();
        }

        /*!\brief Construct as copy from rvalue.*/
        inline bounded(bounded<Type> &&b) noexcept
            : boundedType_(std::move(b.boundedType_)),
              isLowerBounded_(std::move(b.isLowerBounded_)),
              isUpperBounded_(std::move(b.isUpperBounded_)),
              lowerBoundedPtr_(std::move(b.lowerBoundedPtr_)),
              upperBoundedPtr_(std::move(b.upperBoundedPtr_)) {}

        /**
         * \brief Destructor.
         */
        inline ~bounded() noexcept { clear(); }

        inline void clear() noexcept {
            lowerBoundedPtr_.clear();
            upperBoundedPtr_.clear();
        }

        /**
         * \brief Lower bounded.
         */
        hur_nodiscard static inline bounded<Type> lowerBound(const Type &lbv) {
            bounded<Type> lb;
            lb.isLowerBounded_ = true;
            lb.isUpperBounded_ = false;
            lb.lowerBoundedPtr_.reset(new Type(lbv));
            return lb;
        }

        /**
         * \brief Lower bounded.
         */
        hur_nodiscard static inline bounded<Type> lowerBound(const Type &lbv,
                                                             const boundedTypes boundedType) {
            bounded<Type> lb;
            lb.boundedType_ = boundedType;
            lb.isLowerBounded_ = true;
            lb.isUpperBounded_ = false;
            lb.lowerBoundedPtr_.reset(new Type(lbv));
            return lb;
        }

        /**
         * \brief Upper bounded.
         */
        hur_nodiscard static inline bounded<Type> upperBound(const Type &ubv) {
            bounded<Type> lb;
            lb.isLowerBounded_ = false;
            lb.isUpperBounded_ = true;
            lb.upperBoundedPtr_.reset(new Type(ubv));
            return lb;
        }

        /**
         * \brief Upper bounded.
         */
        hur_nodiscard static inline bounded<Type> upperBound(const Type &ubv,
                                                             const boundedTypes boundedType) {
            bounded<Type> lb;
            lb.boundedType_ = boundedType;
            lb.isLowerBounded_ = false;
            lb.isUpperBounded_ = true;
            lb.upperBoundedPtr_.reset(new Type(ubv));
            return lb;
        }

        /**
         * \brief Lower and upper bounded.
         */
        hur_nodiscard static inline bounded<Type>
        lowerUpperBound(const Type &lbv, const Type &ubv,
                        const boundedTypes boundedType = BOUNDED_WITH_VALUE) {
            bounded<Type> lb;
            lb.boundedType_ = boundedType;
            lb.isLowerBounded_ = true;
            lb.isUpperBounded_ = true;
            lb.lowerBoundedPtr_.reset(new Type(lbv));
            lb.upperBoundedPtr_.reset(new Type(ubv));
            return lb;
        }

        /** \brief Is data lower bounded? */
        hur_nodiscard inline bool isLowerBounded() const noexcept { return isLowerBounded_; }

        /** \brief Is data upper bounded? */
        hur_nodiscard inline bool isUpperBounded() const noexcept { return isUpperBounded_; }

        /** \brief The lower bounded value. */
        hur_nodiscard inline const Type &lowerBounded() const noexcept {
            if (!lowerBoundedPtr_) {
                LFatal("Attempt to get access with null lower bounded pointer");
            }
            return *lowerBoundedPtr_;
        }

        /** \brief The upper bounded value. */
        hur_nodiscard inline const Type &upperBounded() const noexcept {
            if (!upperBoundedPtr_) {
                LFatal("Attempt to get access with null lower bounded pointer");
            }
            return *upperBoundedPtr_;
        }

        hur_nodiscard inline boundedTypes boundType() const noexcept { return boundedType_; }

        void checkBoundValue() const {}

        hur_nodiscard inline Type bounding(const Type &v, const Type &limit) const {
            LFatal("Function not define");
        }
    };
} // namespace OpenHurricane
