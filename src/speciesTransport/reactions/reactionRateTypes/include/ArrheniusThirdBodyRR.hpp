/*!
 * \file ArrheniusThirdBodyRR.hpp
 * \brief Header of third-body Arrhenius reaction rate.
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
#include "thirdBodyEfficiency.hpp"
namespace OpenHurricane {
    /**
     * \brief The class of third-body Arrhenius reaction rate.
     */
    class ArrheniusThirdBodyRR : public ArrheniusRR {
    private:
        /**\brief third body efficiencies.*/
        thirdBodyEfficiency thirdBodyEff_;

        /**\brief Concentration of third body*/
        mutable real M_;

    public:
        declareClassName(ArrheniusThirdBodyRR);
        ArrheniusThirdBodyRR() = default;

        inline ArrheniusThirdBodyRR(const ArrheniusRR &ArrCoef, const speciesList &species,
                                    const realArray &efficiencies)
            : ArrheniusRR(ArrCoef), thirdBodyEff_(species, efficiencies), M_(0.0) {}

        inline ArrheniusThirdBodyRR(const real A, const real beta, const real Ta,
                                    const speciesList &species, const realArray &efficiencies)
            : ArrheniusRR(A, beta, Ta), thirdBodyEff_(species, efficiencies), M_(0.0) {}

        inline ArrheniusThirdBodyRR(const ArrheniusThirdBodyRR &atb)
            : ArrheniusRR(atb), thirdBodyEff_(atb.thirdBodyEff_), M_(atb.M_) {}

        ArrheniusThirdBodyRR(const speciesList &sp, const controller &cont);

        virtual inline ~ArrheniusThirdBodyRR() noexcept {}

        hur_nodiscard virtual inline uniquePtr<reactionRateTypes> clone() const {
            return uniquePtr<reactionRateTypes>(new ArrheniusThirdBodyRR(*this));
        }

        hur_nodiscard virtual inline std::string type() const noexcept {
            return std::string("ArrheniusThirdBody");
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
        hur_nodiscard virtual inline real
        DkDT(const real kj, const real p, const real T,
             const realArray &c) const { // The third-body term has been included in k
            return ArrheniusRR::DkDT(kj, p, T, c);
        }

        /**
         * \brief The partial derivatives of third-body terms with respect to the molar concentration of species.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        virtual void gamDGamDCi(const real P, const real T, const realArray &c,
                                realArray &gdgdci) const;

        hur_nodiscard inline virtual bool isModefiedByThirdBody() const noexcept { return true; }

        hur_nodiscard inline virtual bool isArrhenius() const noexcept { return false; }
        hur_nodiscard inline virtual bool isArrheniusThirdBody() const noexcept { return true; }

        hur_nodiscard inline const thirdBodyEfficiency &thirdBodyEff() const noexcept {
            return thirdBodyEff_;
        }
    };
} // namespace OpenHurricane