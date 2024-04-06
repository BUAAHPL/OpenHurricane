/*!
 * \file reactionRateTypes.hpp
 * \brief Header of base class of reaction rate.
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
#include "objectFactory.hpp"
#include "realArray.hpp"
#include "smartPointer.hpp"
#include "speciesList.hpp"

namespace OpenHurricane {
    /**
     * \brief The base class of reaction rate.
     */
    class reactionRateTypes {
    private:
    public:
        declareClassName(reactionRateTypes);
        declareObjFty(reactionRateTypes, controller,
                      (const speciesList &sp, const controller &cont), (sp, cont));

        reactionRateTypes() = default;

        inline reactionRateTypes(const reactionRateTypes &other) {}

        inline reactionRateTypes(const speciesList &sp, const controller &cont) {}

        hur_nodiscard static uniquePtr<reactionRateTypes> creator(const speciesList &sp,
                                                                  const controller &cont);

        virtual ~reactionRateTypes() noexcept {}

        hur_nodiscard virtual inline uniquePtr<reactionRateTypes> clone() const {
            return uniquePtr<reactionRateTypes>(nullptr);
        }

        hur_nodiscard virtual inline std::string type() const noexcept = 0;

        /**
         * \brief Calculating reaction rate.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         * \return The reaction rate conatants.
         */
        hur_nodiscard virtual real k(const real p, const real T, const realArray &c) const = 0;

        /**
         * \brief The partial derivatives of k with respect to temperature.
         * \param[in] kj - The reaction rate conatants.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        hur_nodiscard virtual inline real DkDT(const real kj, const real p, const real T,
                                               const realArray &c) const {
            return 0;
        }

        /**
         * \brief The partial derivatives of k with respect to temperature.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        hur_nodiscard virtual inline real DkDT(const real p, const real T,
                                               const realArray &c) const {
            return DkDT(k(p, T, c), p, T, c);
        }

        /**
         * \brief The partial derivatives of k with respect to pressure.
         * \param[in] kj - The reaction rate conatants.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        hur_nodiscard virtual inline real DkDP(const real kj, const real p, const real T,
                                               const realArray &c) const {
            return real(0);
        }
        /**
         * \brief The partial derivatives of k with respect to pressure.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        hur_nodiscard virtual inline real DkDP(const real p, const real T,
                                               const realArray &c) const {
            return DkDP(k(p, T, c), p, T, c);
        }

        /**
         * \brief The partial derivatives of third-body terms with respect to the molar concentration of species.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        virtual inline void gamDGamDCi(const real P, const real T, const realArray &c,
                                       realArray &gdgdci) const {
            gdgdci.setZero();
        }

        /**
         * \brief The partial derivatives of third-body terms with respect to the temperature.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        hur_nodiscard virtual inline real gamDGamDT(const real P, const real T,
                                                    const realArray &c) const {
            return 0;
        }

        hur_nodiscard inline virtual bool isModefiedByThirdBody() const noexcept { return false; }
        hur_nodiscard inline virtual bool isPressureDenpendent() const noexcept { return false; }
        hur_nodiscard inline virtual bool isArrhenius() const noexcept { return false; }

        hur_nodiscard inline virtual bool isArrheniusThirdBody() const noexcept { return false; }

        hur_nodiscard inline virtual bool isChemicallyActicatedBimolecular() const noexcept {
            return false;
        }

        hur_nodiscard inline virtual bool isUnimolecularFallOff() const noexcept { return false; }
    };
} // namespace OpenHurricane