/*!
 * \file collisionRR.hpp
 * \brief Header of case Collision Frequency Efficiency Expression Rate.
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
#include "ArrheniusRR.hpp"

namespace OpenHurricane {

    /*!\brief The class of bi-molecular spherical collisions frequency efficiency.*/
    class collisionFrequencies {
    private:
        // Private data

        /*!\brief The index of first molecular in specie table.*/
        integer index1_;

        /*!\brief The index of second molecular in specie table.*/
        integer index2_;

        /*!\brief The const reference to the specie table.*/
        const speciesList &species_;

        /*!\brief The average diameter of the two spherical particles.*/
        real d_;

        /*!\brief The reduced molar mass of specie A and B.*/
        real WAB_;

        /*!
         * \brief The reduced gas constant of specie A and B.
         *  i.e. RAB_ = Ru/WAB_.
         */
        real RAB_;

        // Private Member Functions

        void calcAB();

    public:
        inline collisionFrequencies(const integer index1, const integer index2,
                                    const speciesList &sp, const real d)
            : index1_(index1), index2_(index2), species_(sp), d_(d) {
            calcAB();
        }

        inline collisionFrequencies(const collisionFrequencies &other)
            : index1_(other.index1_), index2_(other.index2_), species_(other.species_),
              d_(other.d_), WAB_(other.WAB_), RAB_(other.RAB_) {}

        inline collisionFrequencies(const speciesList &sp, const controller &cont)
            : index1_(cont.findType<integer>("index1", integer())),
              index2_(cont.findType<integer>("index2", integer())), species_(sp),
              d_(cont.findType<real>("diameter", real())) {
            calcAB();
        }

        inline ~collisionFrequencies() noexcept {}

        /*!\brief Return const access to the specie table.*/
        hur_nodiscard inline const speciesList &species() const { return species_; }

        /*!\brief Return the reduced molar mass of specie A and B.*/
        hur_nodiscard inline real WAB() const { return WAB_; }

        /*!\brief Return the reduced gas constant of specie A and B.*/
        hur_nodiscard inline real RAB() const { return RAB_; }

        /*!\brief The basic collision rate.*/
        hur_nodiscard inline real Zb(const real T) const {
            return OpenHurricane::constant::physicalConstant::NA * sqr(d_) *
                   std::sqrt(real(8) * OpenHurricane::constant::mathConstants::pi * RAB_ * T);
        }

        /*!\brief The bi-molecular collision rate.*/
        hur_nodiscard inline real ZAB(const real T, const realArray &c) const {
            return c[index1_] * c[index2_] * Zb(T);
        }

        hur_nodiscard inline real dZbdT(const real T) const {
            return OpenHurricane::constant::physicalConstant::NA * sqr(d_) *
                   std::sqrt(real(2.0) * OpenHurricane::constant::mathConstants::pi * RAB_ / T);
        }
    };

    /**
     * \brief The class of Collision Frequency Efficiency Expression.
     */
    class collisionRR : public reactionRateTypes {
    private:
        collisionFrequencies collFre_;

        ArrheniusRR k_;

    public:
        declareClassName(collisionRR);
        inline collisionRR(const integer index1, const integer index2, const speciesList &sp,
                           const real d, const real A, const real beta, const real Ta)
            : reactionRateTypes(), collFre_(index1, index2, sp, d), k_(A, beta, Ta) {}

        inline collisionRR(const integer index1, const integer index2, const speciesList &sp,
                           const real d, const ArrheniusRR &ka)
            : reactionRateTypes(), collFre_(index1, index2, sp, d), k_(ka) {}

        inline collisionRR(const collisionRR &other)
            : reactionRateTypes(other), collFre_(other.collFre_), k_(other.k_) {}

        inline collisionRR(const speciesList &sp, const controller &cont)
            : reactionRateTypes(sp, cont), collFre_(sp, cont), k_(sp, cont) {}

        virtual inline ~collisionRR() noexcept {}

        hur_nodiscard virtual inline uniquePtr<reactionRateTypes> clone() const {
            return uniquePtr<reactionRateTypes>(new collisionRR(*this));
        }

        hur_nodiscard virtual inline std::string type() const noexcept {
            return std::string("collision");
        }

        /*!\brief Correction factor.*/
        hur_nodiscard inline real gammai(const real p, const real T, const realArray &c) const {
            return OpenHurricane::min(real(1), k_.k(p, T, c));
        }

        /**
         * \brief Calculating reaction rate.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         * \return The reaction rate conatants.
         */
        hur_nodiscard virtual real k(const real p, const real T, const realArray &c) const;

        /**
         * \brief The partial derivatives of k with respect to temperature.
         * \param[in] kj - The reaction rate conatants.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        hur_nodiscard virtual inline real DkDT(const real kj, const real p, const real T,
                                               const realArray &c) const {
            return (kj > 1.0 ? k_.DkDT(kj, p, T, c) : 0.0) * collFre_.dZbdT(T);
        }

        hur_nodiscard inline ArrheniusRR &kArr() noexcept { return k_; }
        hur_nodiscard inline const ArrheniusRR &kArr() const noexcept { return k_; }
    };
} // namespace OpenHurricane