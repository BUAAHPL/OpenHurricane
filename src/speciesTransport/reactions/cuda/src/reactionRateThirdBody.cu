/*!
 * \file reactionRateThirdBody.cu
 * \brief The subroutines and functions of reaction third-body in CUDA platform.
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

#include "reactionRateThirdBody.hpp"
#include <cmath>
#ifdef CUDA_PARALLEL
cu_device cu_real
OpenHurricane::cuChem::thirdBodyCoefficients::thirdBodyEfficiencies(
    const cu_real *__restrict__ c, const cu_ushort irc,
    const cu_ushort nsp) const {
    const auto rii = index_(irc);
    if (rii == -1) {
        return cu_real(1);
    }
    cu_real m = cu_real(0);
    for (cu_ushort isp = 0; isp < nsp; ++isp) {
        m += coefThird_(isp, rii) * c[isp];
    }
    return m;
}

#endif // CUDA_PARALLEL