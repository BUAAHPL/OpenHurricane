/*!
 * \file vectors.hpp
 * \brief Header of vectors.
 *       The subroutines and functions are in the <i>vectors.cpp</i> file.
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

#include "basicFunctions.hpp"
#include "real.hpp"
#include "string.hpp"
#include "vector2DTmpl.hpp"
#include "vectorTmpl.hpp"
#include "complex.hpp"

namespace OpenHurricane {

    using vector = Vector<real>;
    using vector2D = Vector2D<real>;

    template <> inline std::string toString(const vector &is) {
        std::stringstream sx;
        std::stringstream sy;
        std::stringstream sz;
        sx << std::setprecision(feature<real>::precision) << is.x();
        sy << std::setprecision(feature<real>::precision) << is.y();
        sz << std::setprecision(feature<real>::precision) << is.z();
        std::string s;
        s = "(" + sx.str() + ", ";
        s += sy.str() + ", ";
        s += sz.str() + ")";

        return s;
    }

    template <class Type> class flux : public innerProduct<vector, Type> {};

    template <> class flux<real> {
    public:
        typedef real type;
    };

    template <> inline std::string toString(const vector2D &is) {
        std::string s;
        s = "(" + std::to_string(is.x()) + ", ";
        s += std::to_string(is.y()) + ", ";
        return s;
    }
    
    using integerVector = Vector<integer>;
    using integerVector2D = Vector2D<integer>;

    template <> inline std::string toString(const integerVector &is) {
        std::string s;
        s = "(" + std::to_string(is.x()) + ", ";
        s += std::to_string(is.y()) + ", ";
        s += std::to_string(is.z()) + ")";
        return s;
    }

    template <> inline std::string toString(const integerVector2D &is) {
        std::string s;
        s = "(" + std::to_string(is.x()) + ", ";
        s += std::to_string(is.y()) + ")";
        return s;
    }

    using floatVector = Vector<floatReal>;
    using floatVector2D = Vector2D<floatReal>;
#ifdef HURRICANE_DP

    template <> inline std::string toString(const floatVector &is) {
        std::string s;
        s = "(" + std::to_string(is.x()) + ", ";
        s += std::to_string(is.y()) + ", ";
        s += std::to_string(is.z()) + ")";
        return s;
    }

    template <> inline std::string toString(const floatVector2D &is) {
        std::string s;
        s = "(" + std::to_string(is.x()) + ", ";
        s += std::to_string(is.y()) + ", ";
        return s;
    }
#endif // HURRICANE_SP

    using complexVector = Vector<complex>;
    using complexVector2D = Vector2D<complex>;

    template <> hur_nodiscard inline std::string toString(const complexVector &is) {
        std::string s;
        s = "(" + OpenHurricane::toString(is.x()) + ", ";
        s += OpenHurricane::toString(is.y()) + ", ";
        s += OpenHurricane::toString(is.z()) + ")";
        return s;
    }

    template <> hur_nodiscard inline std::string toString(const complexVector2D &is) {
        std::string s;
        s = "(" + OpenHurricane::toString(is.x()) + ", ";
        s += OpenHurricane::toString(is.y()) + ", ";
        return s;
    }

    hur_nodiscard inline complexVector operator*(const complex &v1, const complexVector &v2) {
        return complexVector(v1 * v2.x(), v1 * v2.y(), v1 * v2.z());
    }

    hur_nodiscard inline complexVector operator*(const complexVector &v2, const complex &v1) {
        return complexVector(v1 * v2.x(), v1 * v2.y(), v1 * v2.z());
    }

    hur_nodiscard inline complexVector operator/(const complexVector &v1, const complex &v2) {
        return complexVector(v1.x() / v2, v1.y() / v2, v1.z() / v2);
    }

    hur_nodiscard inline complexVector operator/(const complex &v1, const complexVector &v2) {
        return complexVector(v1 / v2.x(), v1 / v2.y(), v1 / v2.z());
    }

    // complexVector dot product

    hur_nodiscard inline complex operator&(const complexVector &v1, const complexVector &v2) {
        return complex(v1.x() * conjugate(v2.x()) + v1.y() * conjugate(v2.y()) +
                       v1.z() * conjugate(v2.z()));
    }

    // complexVector cross product

    hur_nodiscard inline complexVector operator^(const complexVector &v1, const complexVector &v2) {
        return complexVector((v1.y() * v2.z() - v1.z() * v2.y()),
                             (v1.z() * v2.x() - v1.x() * v2.z()),
                             (v1.x() * v2.y() - v1.y() * v2.x()));
    }

    // complexVector cross product

    hur_nodiscard inline complexVector operator^(const vector &v1, const complexVector &v2) {
        return complexVector((v1.y() * v2.z() - v1.z() * v2.y()),
                             (v1.z() * v2.x() - v1.x() * v2.z()),
                             (v1.x() * v2.y() - v1.y() * v2.x()));
    }

    hur_nodiscard inline complexVector2D operator*(const complex &v1, const complexVector2D &v2) {
        return complexVector2D(v1 * v2.x(), v1 * v2.y());
    }

    hur_nodiscard inline complexVector2D operator*(const complexVector2D &v2, const complex &v1) {
        return complexVector2D(v1 * v2.x(), v1 * v2.y());
    }

    hur_nodiscard inline complexVector2D operator/(const complexVector2D &v1, const complex &v2) {
        return complexVector2D(v1.x() / v2, v1.y() / v2);
    }

    hur_nodiscard inline complexVector2D operator/(const complex &v1, const complexVector2D &v2) {
        return complexVector2D(v1 / v2.x(), v1 / v2.y());
    }

    // complexVector dot product

    hur_nodiscard inline complex operator&(const complexVector2D &v1, const complexVector2D &v2) {
        return complex(v1.x() * conjugate(v2.x()) + v1.y() * conjugate(v2.y()));
    }

    // complexVector cross product

    hur_nodiscard inline complex operator^(const complexVector2D &v1, const complexVector2D &v2) {
        return (v1.x() * v2.y() - v1.y() * v2.x());
    }

    // complexVector cross product

    hur_nodiscard inline complex operator^(const vector2D &v1, const complexVector2D &v2) {
        return (v1.x() * v2.y() - v1.y() * v2.x());
    }
} // namespace OpenHurricane
