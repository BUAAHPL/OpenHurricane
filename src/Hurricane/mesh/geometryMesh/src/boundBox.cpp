/*!
 * \file boundBox.cpp
 * \brief The subroutines and functions of CFD time advance iteration
 * \author Yang Hongzhen
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
#include "boundBox.hpp"

OpenHurricane::boundBox::boundBox(const vector min, const vector max)
    : min_(min), max_(max) {
    for (integer i = 0; i < feature<vector>::nElements_; i++) {
        if (min[i] >= max[i]) {
            LFatal("Min_ of boundBox is greater than Max_.");
        }
    }
}

OpenHurricane::boundBox::boundBox(const vectorArray &points) : min_(), max_() {
    if (points.size() < 2) {
        LFatal("There must be at least two points to construct a boundBox.");
    }
    min_ = points[0];
    max_ = points[0];

    for (integer i = 0; i < points.size(); i++) {
        if (max_.x() < points[i].x()) {
            max_.x() = points[i].x();
        }
        if (max_.y() < points[i].y()) {
            max_.y() = points[i].y();
        }
        if (max_.z() < points[i].z()) {
            max_.z() = points[i].z();
        }

        if (min_.x() > points[i].x()) {
            min_.x() = points[i].x();
        }
        if (min_.y() > points[i].y()) {
            min_.y() = points[i].y();
        }
        if (min_.z() > points[i].z()) {
            min_.z() = points[i].z();
        }
    }
}
