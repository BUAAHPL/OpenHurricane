/*!
 * \file cuChemInterface.hpp
 * \brief Header of interface of reaction from CPU to CUDA platform.
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

namespace OpenHurricane {
    namespace cuChem {
        /**
         * \brief The class of interface of reaction from CPU to CUDA platform.
         */
        class cuChemInterface {
        public:
            enum reactionTypes : short {
                irreversibleReaction = 0,
                reversibleReaction = 1,
                reversibleReactionWithRevPars = 2
            };

            enum thirdBodyTypes : short {
                noThirdBody = 0,         /**< Without third-body. */
                thirdBody = 1,           /**< Only third-body. */
                UnimolecularFallOff = 2, /**< Unimolecular/recombination fall-off reactions. */
                bimolecularFallOff = 3 /**< Chemically activated bimolecular fall-off reactions. */
            };

        public:
            cu_ushort nsp_;
            cu_ushort nrc_;

            /**
             * \brief The class of species list interface to CUDA platform.
             */
            class speciesTableInterface {
            public:
                /** \brief The molecular weight of species */
                cu_real *__restrict__ WPtr_;

                // Temperature limits of applicability of functions
                cu_real *__restrict__ TCommonPtr_;
                cu_real *__restrict__ TLowPtr_;
                cu_real *__restrict__ THighPtr_;

                // Coefficients of JANAF--The pointer of 2D array with size of (7,nsp) [Molar unit]
                cu_real *__restrict__ highCpCoeffsPtr_;
                cu_real *__restrict__ lowCpCoeffsPtr_;

            public:
                /**
                 * \brief Null constructor.
                 */
                inline speciesTableInterface();

                /**
                 * \brief Disallow copy constructor.
                 */
                speciesTableInterface(const speciesTableInterface &) = delete;

                /**
                 * \brief Destructor.
                 */
                inline ~speciesTableInterface() noexcept;

                /**
                 * \brief Transfer.
                 */
                inline void transfer(speciesTableInterface &st);

                /**
                 * \brief Clear all pointer and set null.
                 */
                inline void clear() noexcept;

                inline void alloc(const cu_ushort nsp);
            };

            /**
             * \brief The class of reaction coefficients interface to CUDA platform.
             */
            class reactionCoeffsInterface {
            public:
                cu_ushort maxSpcSizeInRec_;
                cu_ushort *__restrict__ sizePtr_;
                cu_ushort *__restrict__ reactionCoeffsIndexPtr_;
                cu_real *__restrict__ stoichCoeffPtr_;
                cu_real *__restrict__ orderPtr_;

            public:
                /**
                 * \brief Null constructor.
                 */
                inline reactionCoeffsInterface();

                /**
                 * \brief Disallow copy constructor.
                 */
                reactionCoeffsInterface(const reactionCoeffsInterface &) = delete;

                /**
                 * \brief Destructor.
                 */
                inline ~reactionCoeffsInterface() noexcept;

                /**
                 * \brief Transfer.
                 */
                inline void transfer(reactionCoeffsInterface &rc);

                /**
                 * \brief Clear all pointer and set null.
                 */
                inline void clear() noexcept;

                inline void alloc(const cu_ushort nrc, const cu_ushort maxSpcSizeInReac);
            };

            /**
             * \brief The class of reaction rate coefficients interface to CUDA platform.
             */
            class reactionRateCoeffsInterface {
            public:
                cu_real *__restrict__ aPtr_;
                cu_real *__restrict__ bPtr_;
                cu_real *__restrict__ TaPtr_;

            public:
                /**
                 * \brief Null constructor.
                 */
                inline reactionRateCoeffsInterface();

                /**
                 * \brief Disallow copy constructor.
                 */
                reactionRateCoeffsInterface(const reactionRateCoeffsInterface &) = delete;

                /**
                 * \brief Destructor.
                 */
                inline ~reactionRateCoeffsInterface() noexcept;

                /**
                 * \brief Transfer.
                 */
                inline void transfer(reactionRateCoeffsInterface &rc);

                /**
                 * \brief Clear all pointer and set null.
                 */
                inline void clear() noexcept;

                inline void alloc(const cu_integer nrc);
            };

            /**
             * \brief The class of third-body coefficients interface to CUDA platform.
             */
            class thirdBodyCoeffInterface {
            public:
                cu_ushort nrcThird_;
                cu_short *__restrict__ thidrBodyIndexPtr_;
                cu_real *__restrict__ coefThirdPtr_;

            public:
                /**
                 * \brief Null constructor.
                 */
                inline thirdBodyCoeffInterface();

                /**
                 * \brief Disallow copy constructor.
                 */
                thirdBodyCoeffInterface(const thirdBodyCoeffInterface &) = delete;

                /**
                 * \brief Destructor.
                 */
                inline ~thirdBodyCoeffInterface() noexcept;

                /**
                 * \brief Transfer.
                 */
                inline void transfer(thirdBodyCoeffInterface &rc);

                /**
                 * \brief Clear all pointer and set null.
                 */
                inline void clear() noexcept;

                inline void alloc(const cu_ushort nsp, const cu_ushort nrc,
                                  const cu_ushort nrcThird);
            };

            /**
             * \brief The class of pressure dependent reaction coefficients interface to CUDA platform.
             */
            class pressureDepCoeffInterface {
            public:
                enum fallOffTypes {
                    Lindmann = 0,
                    TroeWithoutLastTerm = 1,
                    TroeWithLastTerm = 2,
                    SRI = 3
                };

            public:
                cu_ushort nrcPD_;
                cu_short *__restrict__ indexPtr_;
                cu_real *__restrict__ aPtr_;
                cu_real *__restrict__ bPtr_;
                cu_real *__restrict__ TaPtr_;
                cu_ushort *__restrict__ fallOffTypePtr_;
                cu_real *__restrict__ fallOffCoeffPtr_;

            public:
                /**
                 * \brief Null constructor.
                 */
                inline pressureDepCoeffInterface();

                /**
                 * \brief Disallow copy constructor.
                 */
                pressureDepCoeffInterface(const pressureDepCoeffInterface &) = delete;

                /**
                 * \brief Destructor.
                 */
                inline ~pressureDepCoeffInterface() noexcept;

                /**
                 * \brief Transfer.
                 */
                inline void transfer(pressureDepCoeffInterface &rc);

                /**
                 * \brief Clear all pointer and set null.
                 */
                inline void clear() noexcept;

                inline void alloc(const cu_ushort nrc, const cu_ushort nrcPD);
            };

            /** \brief Species table interface. */
            speciesTableInterface sptInt_;

            /** \brief Forward reaction coefficients interface. */
            reactionCoeffsInterface reacInt_;

            /** \brief Reverse reaction coefficients interface. */
            reactionCoeffsInterface revReacInt_;

            /** \brief The third-body coefficients interface. */
            thirdBodyCoeffInterface thirdBInt_;

            /** \brief The pressure-dependent reaction coefficients interface. */
            pressureDepCoeffInterface pressDepInt_;

            /** \brief Forward reaction rate coefficients interface. */
            reactionRateCoeffsInterface kfInt_;

            /** \brief Reverse reaction rate coefficients interface. */
            reactionRateCoeffsInterface rfInt_;

            /*!\brief The interface of summation of stoichiometry coeffs with size of nrc.*/
            cu_real *__restrict__ sumStoiIntPtr_;

            /**
             * \brief Reaction type:
             *   0 - irreversible reaction.
             *   1 - reversible reaction.
             *   2 - non-equilibrium reversible reaction.
             */
            cu_ushort *__restrict__ reactionTypeIntPtr_;

            /**
             * \brief Third-body type:
             *   0 - Without third-body.
             *   1 - Only third-body.
             *   2 - Unimolecular/recombination fall-off reactions.
             *   3 - Chemically activated bimolecular fall-off reactions.
             */
            cu_ushort *__restrict__ thirdBodyTypeIntPtr_;

        public:
            /**
             * \brief Null constructor.
             */
            inline cuChemInterface();

            /**
             * \brief Disallow copy constructor.
             */
            cuChemInterface(const cuChemInterface &) = delete;

            /**
             * \brief Destructor.
             */
            inline ~cuChemInterface() noexcept;

            /**
             * \brief Transfer.
             */
            inline void transfer(cuChemInterface &rc);

            /**
             * \brief Free pointer.
             */
            inline void clear() noexcept;

            /**
             * \brief Alloc only this level of pointer.
             */
            inline void allocOnlyThis(const cu_ushort nsp, const cu_ushort nrc);
        };
    } // namespace CUDAReactions
} // namespace OpenHurricane

#endif // CUDA_PARALLEL

#include "cuChemInterface.inl"
