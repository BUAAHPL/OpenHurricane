/*!
 * \file reactionRatePressureDependent.hpp
 * \brief Header of pressure dependent reaction in CUDA platform.
 *       The subroutines and functions are in the <i>reactionRatePressureDependent.cpp</i> file.
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
 * \brief The class of pressure dependent coefficients.
 */
class pressureDependentCoefficients {
public:
    enum fallOffTypes {
        Lindmann = 0,
        TroeWithoutLastTerm = 1,
        TroeWithLastTerm = 2,
        SRI = 3
    };

private:
    /*!\brief Index of reactions in the pressure dependent coefficients array.
     *        If the reaction is not pressure dependent reaction, then its index is -1.
     */
    cu1DArray<cu_short> index_;

    /** \brief The pre-exponential factor array with size of nrcPD (number of reactions that are pressure dependent.). */
    cu1DArray<cu_real> a_;

    /** \brief The temperature exponent array with size of nrcPD (number of reactions that are pressure dependent.). */
    cu1DArray<cu_real> b_;

    /** \brief The activation devided by universal gas constant array with size of nrcPD (number of reactions that are pressure dependent.). */
    cu1DArray<cu_real> Ta_;

    /*!\brief The type of fall-off function array with size of nrcPD (number of reactions that are pressure dependent.).
     *        0 - Lindmann.
     *        1 - Troe without last term.
     *        2 - Troe with last term.
     *        3 - SRI
     */
    cu1DArray<cu_ushort> fallOffType_;

    /**
     * \brief The coefficients of fall-off function array with size of nrcPD (number of reactions that are pressure dependent.).
     *        The size is 5 x nrcPD. fallOffCoeff(5,nrcPD)
     * Coefficients:
     *
     * Type|Coeffs
     * -----|------
     * Lindmann|none
     * Troe    |(0,irc)=Alpha,...
     * SRI     |(0,irc)=a,...
     */
    cu2DArray<cu_real> fallOffCoeff_;

public:
    class fallOffFunctions {
    public:
        /**
         * \brief Lindemann fall-off function.
         */
        constexpr inline static cu_dual cu_real
        Lindemann(const cu_real T, const cu_real Pr);

        /**
         * \brief SRI fall-off function.
         */
        static cu_dual cu_real
        SRI(const cu_real a, const cu_real b, const cu_real c,
            const cu_real d, const cu_real e, const cu_real T,
            const cu_real Pr);

        /**
         * \brief Troe fall-off function.
         */
        static cu_dual cu_real Troe(const cu_real alpha,
                                             const cu_real Tsss,
                                             const cu_real Ts,
                                             const cu_real T,
                                             const cu_real Pr);

        /**
         * \brief Troe fall-off function.
         */
        static cu_dual cu_real
        Troe(const cu_real alpha, const cu_real Tsss, const cu_real Ts,
             const cu_real Tss, const cu_real T, const cu_real Pr);
    };

public:
    /**
     * \brief Null constructor.
     */
    inline cu_host pressureDependentCoefficients();

    /**
     * \brief Construct from components.
     * \param[in] nrc - Number of reaction
     * \param[in] nrcPD - Number of reaction that are pressure dependent
     * \param[in] indexPtr - The pointer of index of reactions in the third-body coefficients array with size of nrc
     * \param[in] aPtr - The pointer of pre-exponential factor array with size of nrcPD
     * \param[in] bPtr - The pointer of temperature exponent array with size of nrcPD
     * \param[in] TaPtr - The pointer of activation devided by universal gas constant array with size of nrcPD
     * \param[in] fallOffTypePtr - The pointer of type of fall-off function array with size of nrcPD
     * \param[in] fallOffCoeffPtr - The pointer of coefficients of fall-off function array with size of nrcPD
     */
    inline cu_host pressureDependentCoefficients(
        const cu_ushort nrc, const cu_ushort nrcPD,
        const cu_short *__restrict__ indexPtr,
        const cu_real *__restrict__ aPtr, const cu_real *__restrict__ bPtr,
        const cu_real *__restrict__ TaPtr,
        const cu_ushort *__restrict__ fallOffTypePtr,
        const cu_real *__restrict__ fallOffCoeffPtr);

    /**
     * \brief Construct from components.
     * \param[in] nrc - Number of reaction
     * \param[in] nrcPD - Number of reaction that are pressure dependent
     * \param[in] indexPtr - The pointer of index of reactions in the third-body coefficients array with size of nrc
     * \param[in] aPtr - The pointer of pre-exponential factor array with size of nrcPD
     * \param[in] bPtr - The pointer of temperature exponent array with size of nrcPD
     * \param[in] TaPtr - The pointer of activation devided by universal gas constant array with size of nrcPD
     * \param[in] fallOffTypePtr - The pointer of type of fall-off function array with size of nrcPD
     * \param[in] fallOffCoeffPtr - The pointer of coefficients of fall-off function array with size of nrcPD
     */
    inline cu_host pressureDependentCoefficients(
        const cu_ushort nrc, const cu_ushort nrcPD,
        const cu_short *__restrict__ indexPtr,
        const cu_real *__restrict__ aPtr, const cu_real *__restrict__ bPtr,
        const cu_real *__restrict__ TaPtr,
        const cu_ushort *__restrict__ fallOffTypePtr,
        const cu_real *__restrict__ fallOffCoeffPtr,
        const cudaStreams &streams);

    inline cu_host pressureDependentCoefficients(
        const cu_ushort nrc,
        const cuChemInterface::pressureDepCoeffInterface &pressInt);

    inline cu_host pressureDependentCoefficients(
        const cu_ushort nrc,
        const cuChemInterface::pressureDepCoeffInterface &pressInt,
        const cudaStreams &streams);

    /**
     * \brief Copy constructor.
     */
    inline cu_dual
    pressureDependentCoefficients(const pressureDependentCoefficients &);

    /**
     * \brief Destructor.
     */
    cu_dual inline ~pressureDependentCoefficients() noexcept;

    /**
     * \brief Destroy cuda array memory.
     */
    inline cu_host void destroy();

    /**
     * \brief Return the low pressure limit reaction rate constant for unimolecular/recombination fall-off reactions.
     *        Or high pressure limit reaction rate constant for chemically activated bimolecular fall-off reactions.
     * \param[in] ri - The index of reaction.
     * \param[in] T - Temperature.
     * \return The reaction rate constant.
     * \retval A real value.
     */
    inline cu_device cu_real k(const cu_ushort ri,
                                        const cu_real T) const;

    /**
     * \brief Return the reduced pressure for pressure-dependent reactions.
     * \param[in] k0 - Low pressure limit reaction rate constant.
     * \param[in] kinf - High pressure limit reaction rate constant.
     * \param[in] M - Third-body efficency.
     * \return The reduced pressure.
     * \retval A real value.
     */
    inline cu_device cu_real Pr(const cu_real k0,
                                         const cu_real kinf,
                                         const cu_real M) const;

    /**
     * \brief Return the fall-off function.
     * \param[in] ri - The index of reaction.
     * \param[in] T - Temperature.
     * \param[in] Pr - Reduced pressure.
     * \return The fall-off function.
     * \retval A real value.
     */
    cu_device cu_real F(const cu_ushort ri, const cu_real T,
                                 const cu_real Pr) const;

    /**
     * \brief Unimolecular / Recombination Fall-off Reaction rate.
     * \param[in] kinf - High pressure limit reaction rate constant.
     * \param[in] Pr - The reduced pressure.
     * \param[in] F - Fall-off function.
     * \return The Unimolecular / Recombination Fall-off Reaction rate.
     * \retval A real value.
     */
    inline cu_dual cu_real UnimolecularFallOffRate(
        const cu_real kinf, const cu_real Pr, const cu_real F) const;

    /**
     * \brief Unimolecular / Recombination Fall-off Reaction rate modification factor.
     * \param[in] Pr - The reduced pressure.
     * \param[in] F - Fall-off function.
     * \return The Unimolecular / Recombination Fall-off Reaction rate modification factor.
     * \retval A real value.
     */
    inline cu_dual cu_real
    UnimolecularFallOffFactor(const cu_real Pr, const cu_real F) const;

    /**
     * \brief Chemically Activated Bimolecular Reaction rate.
     * \param[in] k0 - Low pressure limit reaction rate constant.
     * \param[in] Pr - The reduced pressure.
     * \param[in] F - Fall-off function.
     * \return The Chemically Activated Bimolecular Reaction rate.
     * \retval A real value.
     */
    inline cu_dual cu_real BimolecularFallOffRate(
        const cu_real k0, const cu_real Pr, const cu_real F) const;

    /**
     * \brief Chemically Activated Bimolecular Reaction rate modification factor.
     * \param[in] Pr - The reduced pressure.
     * \param[in] F - Fall-off function.
     * \return The Chemically Activated Bimolecular Reaction rate modification factor.
     * \retval A real value.
     */
    inline cu_dual cu_real
    BimolecularFallOffFactor(const cu_real Pr, const cu_real F) const;
};

} // namespace CUDAReactions
} // namespace OpenHurricane

#endif // CUDA_PARALLEL

#include "reactionRatePressureDependent.inl"
