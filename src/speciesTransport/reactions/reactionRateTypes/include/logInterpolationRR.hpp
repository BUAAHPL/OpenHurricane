/*!
 * \file logInterpolationRR.hpp
 * \brief Header of case general pressure dependence reaction rate using logarithmic interpolation.
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
#include "ArrheniusRR.hpp"
#include "thirdBodyEfficiency.hpp"

namespace OpenHurricane {

    /**
     * \brief The class of general pressure dependence reaction rate using logarithmic interpolation.
     */
    class logInterpolationRR : public reactionRateTypes {
    private:
        // Private data

        /**
         * \brief Pressure (unit: [Pa]).
         */
        realArray p_;

        List<ArrheniusRR> k_;

        sharedPtr<thirdBodyEfficiency> thirdBodyEfficiencyPtr_;

        void checkSizeAndPresureRank();

    public:
        declareClassName(logInterpolationRR);

        logInterpolationRR() = default;

        /*!\brief Construct from components without third-body.*/
        inline logInterpolationRR(const realArray &p, const List<ArrheniusRR> &k)
            : reactionRateTypes(), p_(p), k_(k), thirdBodyEfficiencyPtr_(nullptr) {
            checkSizeAndPresureRank();
        }

        /*!\brief Construct from components with third-body.*/
        inline logInterpolationRR(const realArray &p, const List<ArrheniusRR> &k,
                                  const thirdBodyEfficiency &tbe)
            : reactionRateTypes(), p_(p), k_(k),
              thirdBodyEfficiencyPtr_(new thirdBodyEfficiency(tbe)) {
            checkSizeAndPresureRank();
        }

        inline logInterpolationRR(const logInterpolationRR &other)
            : reactionRateTypes(other), p_(other.p_), k_(other.k_),
              thirdBodyEfficiencyPtr_(other.thirdBodyEfficiencyPtr_) {}

        logInterpolationRR(const speciesList &sp, const controller &cont);

        virtual inline ~logInterpolationRR() noexcept {}

        hur_nodiscard virtual inline uniquePtr<reactionRateTypes> clone() const {
            return uniquePtr<reactionRateTypes>(new logInterpolationRR(*this));
        }

        hur_nodiscard virtual inline std::string type() const noexcept {
            return std::string("pressureDependenceUsingLogarithmicInterpolation");
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
         * \brief The partial derivatives of k with respect to pressure.
         * \param[in] kj - The reaction rate conatants.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        hur_nodiscard virtual real DkDP(const real kj, const real p, const real T,
                                        const realArray &c) const;

        /**
         * \brief The partial derivatives of third-body terms with respect to the molar concentration of species.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        virtual void gamDGamDCi(const real P, const real T, const realArray &c,
                                realArray &gdgdci) const;

        hur_nodiscard inline virtual bool isModefiedByThirdBody() const noexcept {
            return !thirdBodyEfficiencyPtr_.isNull();
        }
        hur_nodiscard inline virtual bool isPressureDenpendent() const noexcept { return true; }

        hur_nodiscard inline realArray &ps() noexcept { return p_; }
        hur_nodiscard inline const realArray &ps() const noexcept { return p_; }

        hur_nodiscard inline List<ArrheniusRR> &ks() noexcept { return k_; }
        hur_nodiscard inline const List<ArrheniusRR> &ks() const noexcept { return k_; }

        void resetThirdBodyEff(const thirdBodyEfficiency &tbe);

        inline void append(const real p, const ArrheniusRR &k) {
            p_.append(p);
            k_.append(k);
            checkSizeAndPresureRank();
        }
    };

} // namespace OpenHurricane