/*!
 * \file reactionList.hpp
 * \brief Header of reaction list.
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

#include "reactionTypes.hpp"
#include "thermoList.hpp"

namespace OpenHurricane {

    /**
     * \brief The class of reaction table.
     */
    class reactionList : public sharedPtrList<reactionTypes> {
    public:
        using Base = sharedPtrList<reactionTypes>;

    private:
        // Private data

        /*!\brief Const reference to the species list.*/
        const speciesList &species_;

        /*!\brief Const reference to the thermo table.*/
        const thermoList &thermo_;

        /*!\brief the Gibbs free energy of species at standard pressure.*/
        realArray G0_;

        /*!\brief the Gibbs free energy of species at standard pressure.*/
        realArray DG0DT_;

    public:
        inline reactionList(const speciesList &sp, const thermoList &ther)
            : Base(), species_(sp), thermo_(ther), G0_(sp.size(), Zero), DG0DT_(sp.size(), Zero) {}

        inline reactionList(const speciesList &sp, const thermoList &ther, const controller &cont)
            : Base(), species_(sp), thermo_(ther), G0_(sp.size(), Zero), DG0DT_(sp.size(), Zero) {
            readReactions(cont);
        }

        inline reactionList(const reactionList &other)
            : Base(other), species_(other.species_), thermo_(other.thermo_), G0_(other.G0_),
              DG0DT_(other.DG0DT_) {}

        inline reactionList(reactionList &&other) noexcept
            : Base(std::move(other)), species_(other.species_), thermo_(other.thermo_),
              G0_(std::move(other.G0_)), DG0DT_(std::move(other.DG0DT_)) {}

        virtual ~reactionList() noexcept {}

        void readReactions(const controller &cont);

        /*!
         *\brief Update the standard-state Gibbs free energy field.
         *\param[in] Td - Temperature with dimension K.
         */
        inline void updateG0(const real Td) {
            thermo_.G0(Td, G0_);
            thermo_.DG0DT(Td, DG0DT_);
        }

        hur_nodiscard inline const realArray &G0() const noexcept { return G0_; }

        hur_nodiscard inline realArray &G0() noexcept { return G0_; }

        hur_nodiscard inline const realArray &DG0DT() const noexcept { return DG0DT_; }

        hur_nodiscard inline realArray &DG0DT() noexcept { return DG0DT_; }

        reactionList &operator=(const reactionList &) = delete;
    };

} // namespace OpenHurricane