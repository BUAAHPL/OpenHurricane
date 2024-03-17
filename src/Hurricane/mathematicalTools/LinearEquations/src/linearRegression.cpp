/*!
 * \file linearRegression.cpp
 * \brief Main subroutines for linear regression.
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

#include "linearRegression.hpp"

hur_nodiscard OpenHurricane::real OpenHurricane::slope(const realArray &yy, const realArray &xx) {
    if (yy.size() != xx.size()) {
        LFatal("The size of yy and xx is not the same");
    }
    if (yy.size() == 0) {
        return Zero;
    }
    real sumY = Zero;
    real sumX = Zero;
    real sumXY = Zero;
    real sumX2 = Zero;
    for (integer i = 0; i < yy.size(); ++i) {
        sumXY += (yy[i] * xx[i]);
        sumY += yy[i];
        sumX += xx[i];
        sumX2 += (sqr(xx[i]));
    }

    real aveX = sumX / real(yy.size());
    real aveY = sumY / real(yy.size());

    real tmp = sumX2 - sqr(aveX) * real(yy.size());
    if (tmp) {
        return (sumXY - aveX * aveY * real(yy.size())) / tmp;
    } else {
        return Zero;
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::slope(const realArray &yy) {
    realArray xx(yy.size(), Zero);
    for (integer i = 0; i < yy.size(); ++i) {
        xx[i] = real(i);
    }
    return slope(yy, xx);
}
