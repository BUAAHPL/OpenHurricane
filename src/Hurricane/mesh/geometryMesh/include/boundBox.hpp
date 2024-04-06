/*!
 * \file boundBox.hpp
 * \brief Header of meshInterpolation.
 *       The subroutines and functions are in the <i>boundBox.cpp</i> file.
 * \author Yang Hongzhen
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
#include "vectorArray.hpp"

namespace OpenHurricane {
    class boundBox {
    private:
        /*!\brief bound box min point.*/
        vector min_;

        /*!\brief bound box max point.*/
        vector max_;

    public:
        inline boundBox() : min_(Zero), max_(Zero) {}

        /*!\brief Construct from boundBox min max.*/
        boundBox(const vector min, const vector max);

        /*!\brief Construct from given points.*/
        boundBox(const vectorArray &points);

        inline ~boundBox() noexcept {};

        hur_nodiscard inline const vector &min() const noexcept { return min_; }

        hur_nodiscard inline const vector &max() const noexcept { return max_; }

        hur_nodiscard inline vector &min() noexcept { return min_; }

        hur_nodiscard inline vector &max() noexcept { return max_; }

        hur_nodiscard inline vector span() const { return (max_ - min_); }

        hur_nodiscard inline real mag() const { return (max_ - min_).magnitude(); }

        hur_nodiscard inline real volume() const {
            vector v(max_ - min_);
            return v.x() * v.y() * v.z();
        }

        hur_nodiscard inline bool contains(const vector point) const {
            return (point.x() >= min_.x() && point.x() <= max_.x() && point.y() >= min_.y() &&
                    point.y() <= max_.y() && point.z() >= min_.z() && point.z() <= max_.z());
        }
    };
} // namespace OpenHurricane