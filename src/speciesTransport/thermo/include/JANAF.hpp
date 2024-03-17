/*!
 * \file JANAF.hpp
 * \brief Header of thermo properties by JANAF.
 *       The subroutines and functions are in the <i>JANAF.inl</i> file.
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
#include "controller.hpp"
#include "real.hpp"
#include "thermo.hpp"

namespace OpenHurricane {

    class JANAF : public thermo {
    public:
        /**
         * \brief Array type for coefficients array.
         */
        using coeffArray = FixedList<real, 7>;

    private:
        /**
         * \brief Temperature ranges for 2 sets of coefficients.
         */
        real TLow_, THigh_, TCommon_;

        /**
         * \brief Coefficients array for upper temperature interval.
         */
        coeffArray highCpCoeffs_;

        /**
         * \brief Coefficients array for lower temperature interval.
         */
        coeffArray lowCpCoeffs_;

        /*!\brief Return the coefficients corresponding to the given temperature.*/
        hur_nodiscard inline const coeffArray &coeffs(const real T) const noexcept;

    public:
        declareClassNames;

        /*!
         * \brief Construct from components.
         * \param[in] Tlow - lowest temperature
         * \param[in] Thigh - highest temperature
         * \param[in] Tcommon - common temperature
         * \param[in] highCoeffs - coefficients array for upper temperature interval
         * \param[in] lowCoeffs - coefficients array for lower temperature interval
         */
        inline JANAF(const equationOfState &st, const integer id, const real Tlow, const real Thigh,
                     const real Tcommon, const coeffArray &highCpCoeffs,
                     const coeffArray &lowCpCoeffs, const bool convertCoeffs = false);

        /*!
         * \brief Construct from components.
         * \param[in] Tlow - lowest temperature
         * \param[in] Thigh - highest temperature
         * \param[in] Tcommon - common temperature
         * \param[in] highCoeffs - coefficients array for upper temperature interval
         * \param[in] lowCoeffs - coefficients array for lower temperature interval
         */
        JANAF(const equationOfState &st, const integer id, const real Tlow, const real Thigh,
              const real Tcommon, const coeffArray &highCpCoeffs, const coeffArray &lowCpCoeffs,
              const phaseType pt, const bool convertCoeffs = false);

        /*!
         * \brief Construct from components.
         * \param[in] Tlow - lowest temperature
         * \param[in] Thigh - highest temperature
         * \param[in] Tcommon - common temperature
         * \param[in] highCoeffs - coefficients array for upper temperature interval
         * \param[in] lowCoeffs - coefficients array for lower temperature interval
         */
        JANAF(const equationOfState &st, const integer id, const real Tlow, const real Thigh,
              const real Tcommon, const coeffArray &highCpCoeffs, const coeffArray &lowCpCoeffs,
              const std::string &pt, const bool convertCoeffs = false);

        JANAF(const controller &cont, const equationOfState &st, const integer id);

        JANAF(const JANAF &);
        JANAF &operator=(const JANAF &j);

        JANAF(const JANAF &, const string &);

        /*!\brief Return a clone.*/
        virtual hur_nodiscard hur_nodiscard uniquePtr<thermo> clone() const {
            return uniquePtr<thermo>(new JANAF(*this));
        }

        virtual ~JANAF() noexcept {}

        /*!\brief Limit the temperature to be in the range TLow_ to THigh_.*/
        hur_nodiscard inline real limit(const real T) const noexcept;

        /*!
         * \brief Limit the temperature to be in the range TLow_ to THigh_.
         * \param[out] iFlag -# 0: T is lower than Tlow_;
         *                   -# 1: T is common;
         *                   -# 2: T is greater than THigh_.
         */
        inline real limit(const real T, integer &iFlag) const noexcept;

        /*!\brief Return const access to the low temperature limit.*/
        hur_nodiscard inline real Tlow() const noexcept;

        /*!\brief Return const access to the high temperature limit.*/
        hur_nodiscard inline real Thigh() const noexcept;

        /*!\brief Return const access to the common temperature.*/
        hur_nodiscard inline real Tcommon() const noexcept;

        /*!\brief Return const access to the high temperature poly coefficients.*/
        hur_nodiscard inline const coeffArray &highCpCoeffs() const noexcept;

        /*!\brief Return const access to the low temperature poly coefficients.*/
        hur_nodiscard inline const coeffArray &lowCpCoeffs() const noexcept;

        /*!\brief Heat capacity at constant pressure [J/(kg K)] at standard pressure (1 atm).*/
        virtual hur_nodiscard real cp0(const real T) const noexcept;

        /*!\brief Absolute Enthalpy [J/kg] at standard pressure.*/
        virtual hur_nodiscard real ha0(const real T) const noexcept;

        /*!\brief Absolute Enthalpy [J/kg] at standard pressure.*/
        inline virtual hur_nodiscard real Dha0DT(const real T) const noexcept;

        /*!\brief Chemical enthalpy [J/kg].*/
        virtual hur_nodiscard real hc() const noexcept;

        /*!\brief Entropy [J/(kg K)] at standard pressure.*/
        virtual hur_nodiscard real s0(const real T) const;

        /*!\brief The Gibbs free energy of species at standard pressure.*/
        hur_nodiscard virtual real g0(const real T) const;

        /*!\brief Entropy [J/(kg K)] at standard pressure.*/
        inline virtual hur_nodiscard real Ds0DT(const real T) const noexcept;

        /*!\brief Heat capacity at constant pressure [J/(kg K)].*/
        virtual hur_nodiscard real Dcp_pDT(const real pi, const real T) const noexcept;

    };

} // namespace OpenHurricane

#include "JANAF.inl"