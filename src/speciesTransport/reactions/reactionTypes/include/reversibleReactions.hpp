/*!
 * \file reversibleReactions.hpp
 * \brief Header of reversible reactions.
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

#include "reactionTypes.hpp"
namespace OpenHurricane {

    class reversibleReactions : public reactionTypes {
    private:
        /*!\brief Forward reaction rate.*/
        uniquePtr<reactionRateTypes> kfPtr_;

    public:
        declareClassName(reversibleReactions);

        reversibleReactions(const reactionTypes &reac, const reactionRateTypes &rrf);

        reversibleReactions(const speciesList &sp, const reactionList &rt, const controller &cont);

        reversibleReactions(const reversibleReactions &other);

        hur_nodiscard virtual inline uniquePtr<reactionTypes> clone() const {
            return uniquePtr<reactionTypes>(new reversibleReactions(*this));
        }

        /*!\brief Destructor.*/
        virtual inline ~reversibleReactions() noexcept {}

        // Reaction rate constant

        /*!\brief Forward rate constant.*/
        hur_nodiscard virtual inline real kf(const real p, const real T, const realArray &c) const {
            return kfPtr_->k(p, T, c);
        }

        /*!\brief Reverse rate constant from the given forward constant.*/
        hur_nodiscard virtual inline real kr(const real kfwd, const real p, const real T,
                                             const realArray &c) const {
            const real Kc = this->Kc(T);
            if (mag(Kc) > veryTiny) {
                return kfwd / Kc;
            } else {
                return Zero;
            }
        }

        /*!\brief Reverse rate constant.*/
        hur_nodiscard virtual inline real kr(const real p, const real T, const realArray &c) const {
            return kr(kfPtr_->k(p, T, c), p, T, c);
        }

        hur_nodiscard virtual inline bool hasThirdBodyTerms() const noexcept {
            return kfPtr_->isModefiedByThirdBody();
        }

        hur_nodiscard virtual inline bool isPressureDenpendent() const noexcept {
            return kfPtr_->isPressureDenpendent();
        }

        hur_nodiscard virtual inline bool isModefiedByThirdBody() const noexcept {
            return kfPtr_->isModefiedByThirdBody();
        }

        /*!\brief Return the temperature derivative of forward rate.
         * \param[in] kf - forward rate
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \return Return the temperature derivative of forward rate.
         */
        hur_nodiscard virtual inline real DkfDT(const real kf, const real p, const real T,
                                                const realArray &c) const {
            return kfPtr_->DkDT(kf, p, T, c);
        }

        /*!\brief Return the temperature derivative of forward rate.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \return Return the temperature derivative of forward rate.
         */
        hur_nodiscard virtual inline real DkfDT(const real p, const real T,
                                                const realArray &c) const {
            return kfPtr_->DkDT(p, T, c);
        }

        /*!\brief Return the temperature derivative of forward rate.
         * \param[in] kr - forward rate
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \return Return the temperature derivative of forward rate.
         */
        hur_nodiscard virtual inline real DkrDT(const real kr, const real p, const real T,
                                                const realArray &c, const real dkfdT) const {
            const real Kc = max(this->Kc(T), veryTiny);
            return (dkfdT - kr * DKcDT(T)) / Kc;
        }

        /** \brief The partial derivatives of third-body terms with respect to the species molar concentrations.*/
        virtual void inline DGGDc(const real p, const real T, const realArray &c,
                                  realArray &dgdci) const {
            kfPtr_->gamDGamDCi(p, T, c, dgdci);
        }

        hur_nodiscard virtual inline real DGGDT(const real p, const real T,
                                                const realArray &c) const {
            return kfPtr_->gamDGamDT(p, T, c);
        }

        /*!\brief Return the pressure derivative of forward rate.
         * \param[in] kf - forward rate
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \return Return the temperature derivative of forward rate.
         */
        hur_nodiscard virtual inline real DkfDp(const real kf, const real p, const real T,
                                                const realArray &c) const {
            return kfPtr_->DkDP(kf, p, T, c);
        }

        /*!\brief Return the pressure derivative of forward rate.
         * \param[in] kr - forward rate
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \return Return the temperature derivative of forward rate.
         */
        hur_nodiscard virtual inline real DkrDp(const real kr, const real p, const real T,
                                                const realArray &c, const real dkfdp) const {
            const real Kc = max(this->Kc(T), veryTiny);
            return dkfdp / Kc;
        }
#ifdef CUDA_PARALLEL

    public:
        virtual void cuSetReactionType(const integer reacId, const integer nrc,
                                       cuChem::cuChemInterface &reacInt) const;

#endif // CUDA_PARALLEL
    };

} // namespace OpenHurricane