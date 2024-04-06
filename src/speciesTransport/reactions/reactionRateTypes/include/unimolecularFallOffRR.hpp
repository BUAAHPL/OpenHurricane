/*!
 * \file unimolecularFallOffRR.hpp
 * \brief Header of unimolecular fall-off reaction rate.
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
#include "fallOffFunctions.hpp"
#include "thirdBodyEfficiency.hpp"

namespace OpenHurricane {
    /**
     * \brief The class of unimolecular/recombination fall-off reactions.
     */
    class unimolecularFallOffRR : public reactionRateTypes {
    private:
        // Private data

        /**\brief The low-pressure limit rate.*/
        ArrheniusRR k0_;

        /**\brief The high-pressure limit rate.*/
        ArrheniusRR kinf_;

        /**\brief The blending function.*/
        uniquePtr<fallOffFunctions::fallOffFunction> F_;

        /**\brief Third body efficiencies.*/
        thirdBodyEfficiency thirdBodyEfficiency_;

    public:
        declareClassName(unimolecularFallOffRR);

        inline unimolecularFallOffRR(const ArrheniusRR &k0, const ArrheniusRR &kinf,
                                     const fallOffFunctions::fallOffFunction &F,
                                     const thirdBodyEfficiency &tbe)
            : reactionRateTypes(), k0_(k0), kinf_(kinf), F_(F.clone()), thirdBodyEfficiency_(tbe) {}

        inline unimolecularFallOffRR(const ArrheniusRR &k0, const ArrheniusRR &kinf,
                                     const fallOffFunctions::fallOffFunction &F,
                                     const speciesList &species, const realArray &efficiencies)
            : reactionRateTypes(), k0_(k0), kinf_(kinf), F_(F.clone()),
              thirdBodyEfficiency_(species, efficiencies) {}

        inline unimolecularFallOffRR(const unimolecularFallOffRR &other)
            : reactionRateTypes(other), k0_(other.k0_), kinf_(other.kinf_), F_((*other.F_).clone()),
              thirdBodyEfficiency_(other.thirdBodyEfficiency_) {}

        inline unimolecularFallOffRR(const speciesList &sp, const controller &cont)
            : reactionRateTypes(sp, cont), k0_(sp, cont.subController("low")),
              kinf_(sp, cont.subController("high")), F_(nullptr), thirdBodyEfficiency_(sp, cont) {
            F_ = fallOffFunctions::fallOffFunction::creator(cont.subController("fallOfFunction"));
        }

        inline virtual ~unimolecularFallOffRR() noexcept {}

        hur_nodiscard virtual inline uniquePtr<reactionRateTypes> clone() const {
            return uniquePtr<reactionRateTypes>(new unimolecularFallOffRR(*this));
        }

        hur_nodiscard virtual inline std::string type() const noexcept {
            return std::string("unimolecularFallOff");
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
        hur_nodiscard virtual real DkDT(const real kj, const real p, const real T,
                                        const realArray &c) const;

        /**
         * \brief The partial derivatives of third-body terms with respect to the molar concentration of species.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        virtual void gamDGamDCi(const real P, const real T, const realArray &c,
                                realArray &gdgdci) const;

        /**
         * \brief The partial derivatives of third-body terms with respect to the temperature.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        hur_nodiscard virtual real gamDGamDT(const real P, const real T, const realArray &c) const;

        hur_nodiscard inline virtual bool isModefiedByThirdBody() const noexcept { return true; }

        hur_nodiscard inline virtual bool isUnimolecularFallOff() const noexcept { return true; }

        /**\brief The low-pressure limit rate.*/
        hur_nodiscard inline ArrheniusRR &k0() noexcept { return k0_; }

        /**\brief The low-pressure limit rate.*/
        hur_nodiscard inline const ArrheniusRR &k0() const noexcept { return k0_; }

        /**\brief The high-pressure limit rate.*/
        hur_nodiscard inline ArrheniusRR &kinf() noexcept { return kinf_; }

        /**\brief The high-pressure limit rate.*/
        hur_nodiscard inline const ArrheniusRR &kinf() const noexcept { return kinf_; }

        inline void resetFallOff(const fallOffFunctions::fallOffFunction &fof) {
            F_.clear();
            F_ = fof.clone();
        }

        /**\brief The blending function.*/
        hur_nodiscard inline fallOffFunctions::fallOffFunction &fallOff() noexcept { return *F_; }
        hur_nodiscard inline const fallOffFunctions::fallOffFunction &fallOff() const noexcept {
            return *F_;
        }

        hur_nodiscard inline const thirdBodyEfficiency &thirdBodyEff() const noexcept {
            return thirdBodyEfficiency_;
        }

        void resetThirdBodyEff(const thirdBodyEfficiency &tbe);
    };

} // namespace OpenHurricane
