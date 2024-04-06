/*!
 * \file reactionRateCoeffsCUDA.hpp
 * \brief Header of reaction rate in CUDA platform.
 *       The subroutines and functions are in the <i>reactionRateCoeffsCUDA.cpp</i> file.
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
 * \brief The class of reaction coefficients.
 */
class reactionCoeffs {
public:
    /** \brief The size of species in the coefficients array of each reaction. */
    cu1DArray<cu_ushort> size_;

    /*!\brief Index of species. (maxSize,nrc)*/
    cu2DArray<cu_ushort> index_;

    /*!\brief stoichiometry coefficients. (maxSize,nrc)*/
    cu2DArray<cu_real> stoichCoeff_;

    /*!\brief Reaction order. (maxSize,nrc)*/
    cu2DArray<cu_real> order_;

public:
    /**
     * \brief Null constructor.
     */
    inline cu_host reactionCoeffs();

    /**
     * \brief Construct from components.
     * \param[in] nrc - Number of reactions
     * \param[in] maxSize - The max size of species of each reaction
     * \param[in] sizePtr - The pointer of size array
     * \param[in] indexPtr - The pointer of index array
     * \param[in] stoichCoeffPtr - The pointer of stoichiometry coefficients array
     * \param[in] orderPtr - The pointer of reaction order array
     */
    inline cu_host
    reactionCoeffs(const cu_ushort nrc, const cu_ushort maxSize,
                   const cu_ushort *__restrict__ sizePtr,
                   const cu_ushort *__restrict__ indexPtr,
                   const cu_real *__restrict__ stoichCoeffPtr,
                   const cu_real *__restrict__ orderPtr);

    /**
     * \brief Construct from components.
     * \param[in] nrc - Number of reactions
     * \param[in] maxSize - The max size of species of each reaction
     * \param[in] sizePtr - The pointer of size array
     * \param[in] indexPtr - The pointer of index array
     * \param[in] stoichCoeffPtr - The pointer of stoichiometry coefficients array
     * \param[in] orderPtr - The pointer of reaction order array
     */
    inline cu_host
    reactionCoeffs(const cu_ushort nrc, const cu_ushort maxSize,
                   const cu_ushort *__restrict__ sizePtr,
                   const cu_ushort *__restrict__ indexPtr,
                   const cu_real *__restrict__ stoichCoeffPtr,
                   const cu_real *__restrict__ orderPtr,
                   const cudaStreams &streams);

    inline cu_host reactionCoeffs(
        const cu_ushort nrc,
        const cuChemInterface::reactionCoeffsInterface &reacInt);

    inline cu_host reactionCoeffs(
        const cu_ushort nrc,
        const cuChemInterface::reactionCoeffsInterface &reacInt,
        const cudaStreams &streams);

    /**
     * \brief Copy constructor.
     */
    inline cu_dual reactionCoeffs(const reactionCoeffs &);

    /**
     * \brief Destructor.
     */
    inline cu_dual ~reactionCoeffs() noexcept;

    /**
     * \brief Destroy cuda array memory.
     */
    inline cu_host void destroy();

    /** \brief The size of species in the coefficients array of irc_th reaction. */
    inline cu_device cu_ushort size(const cu_ushort irc) const;
};

/**
 * \brief The class of reaction rate coefficients.
 */
class reactionRateCoeffs {
private:
    /** \brief The pre-exponential factor array with size of nrc. */
    cu1DArray<cu_real> a_;

    /** \brief The temperature exponent array with size of nrc. */
    cu1DArray<cu_real> b_;

    /** \brief The activation devided by universal gas constant array with size of nrc. */
    cu1DArray<cu_real> Ta_;

    /**
     * \brief Return Arrhenius reaction rate in CUDA host and device.
     * \param[in] a - The pre-exponential factor.
     * \param[in] b - The temperature exponent.
     * \param[in] Ta - The activation devided by universal gas constant.
     * \param[in] ri - The index of reaction.
     * \param[in] T - The static temperature.
     * \return The reaction rate conatants.
     * \retval A real value.
     */
    inline cu_dual cu_real ArrheniusReactionRate(
        const cu_real *__restrict__ a, const cu_real *__restrict__ b,
        const cu_real *__restrict__ Ta, const cu_integer ri,
        const cu_real T) const;

    bool destroyed_;

public:
    /**
     * \brief Null constructor.
     */
    inline cu_host reactionRateCoeffs();

    /**
     * \brief Construct from components.
     * \param[in] nrc - Number of reaction
     * \param[in] aPtr - The pointer of pre-exponential factor array
     * \param[in] bPtr - The pointer of temperature exponent array
     * \param[in] TaPtr - The pointer of activation devided by universal gas constant array
     */
    inline cu_host
    reactionRateCoeffs(const cu_ushort nrc,
                       const cu_real *__restrict__ aPtr,
                       const cu_real *__restrict__ bPtr,
                       const cu_real *__restrict__ TaPtr);

    /**
     * \brief Construct from components.
     * \param[in] nrc - Number of reaction
     * \param[in] aPtr - The pointer of pre-exponential factor array
     * \param[in] bPtr - The pointer of temperature exponent array
     * \param[in] TaPtr - The pointer of activation devided by universal gas constant array
     */
    inline cu_host reactionRateCoeffs(
        const cu_ushort nrc, const cu_real *__restrict__ aPtr,
        const cu_real *__restrict__ bPtr, const cu_real *__restrict__ TaPtr,
        const cudaStreams &streams);

    inline cu_host reactionRateCoeffs(
        const cu_ushort nrc,
        const cuChemInterface::reactionRateCoeffsInterface &reacInt);

    inline cu_host reactionRateCoeffs(
        const cu_ushort nrc,
        const cuChemInterface::reactionRateCoeffsInterface &reacInt,
        const cudaStreams &streams);

    /**
     * \brief Copy constructor.
     */
    inline cu_dual reactionRateCoeffs(const reactionRateCoeffs &);

    /**
     * \brief Destructor.
     */
    cu_dual inline ~reactionRateCoeffs() noexcept;

    /**
     * \brief Destroy cuda array memory.
     */
    inline cu_host void destroy();

    inline cu_device cu_real k(const cu_ushort ri,
                                        const cu_real T) const;
};
} // namespace CUDAReactions
} // namespace OpenHurricane

#endif // CUDA_PARALLEL

#include "reactionRateCoeffsCUDA.inl"
