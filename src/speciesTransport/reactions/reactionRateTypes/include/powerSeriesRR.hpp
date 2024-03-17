/*!
 * \file powerSeriesRR.hpp
 * \brief Header of case power series reaction rate.
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
#include "FixedList.hpp"
#include "reactionRateTypes.hpp"

namespace OpenHurricane {
    /**
     * \brief The class of power series reaction rate constants.
     */
    class powerSeriesRR : public reactionRateTypes {
    private:
        // Private data
        real A_;
        real beta_;

        static constexpr integer nb_ = 4;
        FixedList<real, nb_> b_;

    public:
        declareClassName(powerSeriesRR);

        inline powerSeriesRR(const real A, const real beta, const FixedList<real, nb_> &b)
            : reactionRateTypes(), A_(A), beta_(beta), b_(b) {}

        inline powerSeriesRR(const powerSeriesRR &other)
            : reactionRateTypes(other), A_(other.A_), beta_(other.beta_), b_(other.b_) {}
        inline powerSeriesRR(const speciesList &sp, const controller &cont)
            : reactionRateTypes(sp, cont), A_(cont.findType<real>("A", real(0.0))),
              beta_(cont.findType<real>("beta", real(0.0))),
              b_(cont.findType<FixedList<real, nb_>>("b", FixedList<real, nb_>())) {}

        virtual inline ~powerSeriesRR() noexcept {}

        hur_nodiscard virtual inline uniquePtr<reactionRateTypes> clone() const {
            return uniquePtr<reactionRateTypes>(new powerSeriesRR(*this));
        }

        hur_nodiscard virtual inline std::string type() const noexcept {
            return std::string("powerSeries");
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
         * \brief The pre-exponential factor..
         */
        hur_nodiscard inline real &A() noexcept { return A_; }
        hur_nodiscard inline real A() const noexcept { return A_; }

        /**
         * \brief The temperature exponent.
         */
        hur_nodiscard inline real &beta() noexcept { return beta_; }
        hur_nodiscard inline real beta() const noexcept { return beta_; }

        hur_nodiscard inline FixedList<real, nb_> &b() noexcept { return b_; }
        hur_nodiscard inline FixedList<real, nb_> b() const noexcept { return b_; }
    };

} // namespace OpenHurricane