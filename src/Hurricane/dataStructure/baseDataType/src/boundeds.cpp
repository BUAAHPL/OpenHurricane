/*!
 * \file boundeds.cpp
 * \brief Main subroutines for boundeds.
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

#include "boundeds.hpp"

template <> void OpenHurricane::realBounded::checkBoundValue() const {
    if (isLowerBounded_ && isUpperBounded_) {
        if (lowerBounded() > upperBounded()) {
            errorAbortStr(
                ("The lower baounded value: " + toString(lowerBounded()) +
                 " is greater than the upper bounded value: " + toString(upperBounded())));
        }
    }
}

template <> void OpenHurricane::vectorBounded::checkBoundValue() const {
    if (isLowerBounded_ && isUpperBounded_) {
        if (lowerBounded().x() > upperBounded().x() || lowerBounded().y() > upperBounded().y() ||
            lowerBounded().z() > upperBounded().z()) {
            errorAbortStr(
                ("The lower baounded value: " + toString(lowerBounded()) +
                 " is greater than the upper bounded value: " + toString(upperBounded())));
        }
    }
}

template <>
hur_nodiscard OpenHurricane::real
OpenHurricane::bounded<OpenHurricane::real>::bounding(const real &v, const real &limit) const {
    if (isLowerBounded_ && isUpperBounded_) {
        if (boundedType_ == BOUNDED_WITH_VALUE) {
            return min(max(v, *lowerBoundedPtr_), *upperBoundedPtr_);
        } else {
            if (v < *lowerBoundedPtr_ || v > *upperBoundedPtr_) {
                return limit;
            }
        }
    } else if (isLowerBounded_) {
        if (boundedType_ == BOUNDED_WITH_VALUE) {
            return max(v, *lowerBoundedPtr_);
        } else if (v < *lowerBoundedPtr_) {
            return limit;
        }
    } else if (isUpperBounded_) {
        if (boundedType_ == BOUNDED_WITH_VALUE) {
            return min(v, *upperBoundedPtr_);
        } else if (v > *upperBoundedPtr_) {
            return limit;
        }
    }
    return v;
}

template <>
hur_nodiscard OpenHurricane::vector
OpenHurricane::bounded<OpenHurricane::vector>::bounding(const vector &v, const vector &limit) const {
    if (isLowerBounded_ && isUpperBounded_) {
        if (boundedType_ == BOUNDED_WITH_VALUE) {
            return componentMin(componentMax(v, *lowerBoundedPtr_), *upperBoundedPtr_);
        } else {
            if (v.x() < lowerBoundedPtr_->x() || v.x() > upperBoundedPtr_->x() ||
                v.y() < lowerBoundedPtr_->y() || v.y() > upperBoundedPtr_->y() ||
                v.z() < lowerBoundedPtr_->z() || v.z() > upperBoundedPtr_->z()) {
                return limit;
            }
        }
    } else if (isLowerBounded_) {
        if (boundedType_ == BOUNDED_WITH_VALUE) {
            return componentMax(v, *lowerBoundedPtr_);
        } else if (v.x() < lowerBoundedPtr_->x() || v.y() < lowerBoundedPtr_->y() ||
                   v.z() < lowerBoundedPtr_->z()) {
            return limit;
        }
    } else if (isUpperBounded_) {
        if (boundedType_ == BOUNDED_WITH_VALUE) {
            return componentMin(v, *upperBoundedPtr_);
        } else if (v.x() > upperBoundedPtr_->x() || v.y() > upperBoundedPtr_->y() ||
                   v.z() > upperBoundedPtr_->z()) {
            return limit;
        }
    }
    return v;
}
