/*!
 * \file reactionRateThirdBody.hpp
 * \brief Header of reaction rate third-body in CUDA platform.
 *       The subroutines and functions are in the <i>reactionRateThirdBody.cpp</i> file.
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

#ifdef CUDA_PARALLEL
#include "CUDAFunctions.hpp"
#include "JANAFCUDA.hpp"
#include "cudaStreams.hpp"
#include "physicalThermoConatsnts.hpp"
#include "cuChemInterface.hpp"
#include <cmath>

namespace OpenHurricane {
namespace cuChem {
/**
 * \brief The class of third-body coefficients.
 */
class thirdBodyCoefficients {
private:
    /*!\brief Index of reactions in the third-body coefficients array(nrc).
     *        If the reaction is without third-body, then its index is -1.
     */
    cu1DArray<cu_short> index_;

    /*!\brief The third-body coefficients array stored as: coefThird(nsp,nrcThird).
     */
    cu2DArray<cu_real> coefThird_;

public:
    /**
     * \brief Null constructor.
     */
    inline cu_host thirdBodyCoefficients();

    /**
     * \brief Construct from components.
     * \param[in] nsp - Number of species
     * \param[in] nrc - Number of reaction
     * \param[in] nrcThird - Number of reaction
     * \param[in] indexPtr - The pointer of index of reactions in the third-body coefficients array
     * \param[in] coefThirdPtr - The pointer of third-body coefficients array stored as: coefThird(nsp,nrcThird).
     *							 Where, nrcThird is the number of reactions that contain the third-body
     */
    inline cu_host
    thirdBodyCoefficients(const cu_ushort nsp, const cu_ushort nrc,
                          const cu_ushort nrcThird,
                          const cu_short *__restrict__ indexPtr,
                          const cu_real *__restrict__ coefThirdPtr);

    /**
     * \brief Construct from components.
     * \param[in] nsp - Number of species
     * \param[in] nrc - Number of reaction
     * \param[in] nrcThird - Number of reaction
     * \param[in] indexPtr - The pointer of index of reactions in the third-body coefficients array
     * \param[in] coefThirdPtr - The pointer of third-body coefficients array stored as: coefThird(nsp,nrcThird).
     *							 Where, nrcThird is the number of reactions that contain the third-body
     */
    inline cu_host thirdBodyCoefficients(
        const cu_ushort nsp, const cu_ushort nrc,
        const cu_ushort nrcThird, const cu_short *__restrict__ indexPtr,
        const cu_real *__restrict__ coefThirdPtr, const cudaStreams &streams);

    inline cu_host thirdBodyCoefficients(
        const cu_ushort nsp, const cu_ushort nrc,
        const cuChemInterface::thirdBodyCoeffInterface &thirdBInt);

    inline cu_host thirdBodyCoefficients(
        const cu_ushort nsp, const cu_ushort nrc,
        const cuChemInterface::thirdBodyCoeffInterface &thirdBInt,
        const cudaStreams &streams);

    /**
     * \brief Copy constructor.
     */
    inline cu_dual thirdBodyCoefficients(const thirdBodyCoefficients &);

    /**
     * \brief Destructor.
     */
    cu_dual inline ~thirdBodyCoefficients() noexcept;

    /**
     * \brief Destroy cuda array memory.
     */
    inline cu_host void destroy();

    /**
     * \brief Compute the third-body efficiencies.
     * \param[i] c - The pointer of species molar concentrations
     * \param[in] irc - The index of the reaction to calculate the third-body efficiencies
     * \param[in] nsp - The size of species
     * \return The third-body efficiencies.
     * \retval A real value.
     */
    cu_device cu_real
    thirdBodyEfficiencies(const cu_real *__restrict__ c,
                          const cu_ushort irc, const cu_ushort nsp) const;
};
} // namespace CUDAReactions
} // namespace OpenHurricane

#endif // CUDA_PARALLEL

#include "reactionRateThirdBody.inl"
