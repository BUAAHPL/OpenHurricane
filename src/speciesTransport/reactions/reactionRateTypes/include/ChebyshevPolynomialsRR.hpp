/*!
 * \file ChebyshevPolynomialsRR.hpp
 * \brief Header of multiple-well multiple-channel reactions using Chebyshev polynomials.
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
#include "reactionRateTypes.hpp"
#include "thirdBodyEfficiency.hpp"

namespace OpenHurricane {
    /**
     * \brief The namespace of pressure limits and tempetrature limits in Chebyshev polynomials.
     */
    namespace ChebyshevRRLimit {
        /**
         * \brief The namespace of pressure limits in Chebyshev polynomials.
         */
        namespace pressureLimit {
            /**\brief Default minimun pressure limit (Unit: Pa).*/
            static real Pmin = 0.001 * constant::physicalConstant::Patm;

            /**\brief Default maximun pressure limit (Unit: Pa).*/
            static real Pmax = 100.0 * constant::physicalConstant::Patm;
        } // namespace pressureLimit

        /**
         * \brief The namespace of tempetrature limits in Chebyshev polynomials.
         */
        namespace temperatureLimit {
            /**\brief Default minimun temperature limit (Unit: K).*/
            static real Tmin = 300.0;

            /**\brief Default maximun temperature limit (Unit: K).*/
            static real Tmax = 2500.0;
        } // namespace temperatureLimit

    } // namespace ChebyshevRRLimit

    /**
     * \brief The class of Multiple-well Multiple-channel reactions using Chebyshev polynomials.
     */
    class ChebyshevPolynomialsRR : public reactionRateTypes {
    private:
        // Private data

        /**\brief The n x m coefficients.*/
        realArrayArray Anm_;

        real Tmin_, powTminusOneMin_;
        real Tmax_, powTminusOneMax_;

        real Pmin_, log10Pmin_;
        real Pmax_, log10Pmax_;

        uniquePtr<thirdBodyEfficiency> thirdBodyEfficiency_;

        /**
         *\brief Chebyshev polynomials:
         *       n = 0,1,...
         *      Î¦(x) = cos(n*arccos(x)), n = 0,1,...
         */
        hur_nodiscard inline real phi(const integer n, const real x) const {
            return std::cos(real(n) * std::acos(x));
        }

        hur_nodiscard inline real DPhiDX(const integer n, const real x) const {
            return std::sin(real(n - 1) * std::acos(x)) * real(n - 1) /
                   sqrt(max(real(1) - sqr(x), tiny));
        }

        /**\brief The x domain is mapped onto a square bounded by +-1 using this transformtion.*/
        hur_nodiscard inline real transformation(const real x, const real minx,
                                                 const real maxx) const {
            return (2.0 * x - minx - maxx) / (maxx - minx);
        }

        realArray tmpAnm_;

    public:
        declareClassName(ChebyshevPolynomialsRR);

        inline ChebyshevPolynomialsRR()
            : reactionRateTypes(), Anm_(), Tmin_(ChebyshevRRLimit::temperatureLimit::Tmin),
              Tmax_(ChebyshevRRLimit::temperatureLimit::Tmax),
              Pmin_(ChebyshevRRLimit::pressureLimit::Pmin),
              Pmax_(ChebyshevRRLimit::pressureLimit::Pmax), thirdBodyEfficiency_(nullptr) {
            powTminusOneMin_ = 1.0 / Tmin_;
            powTminusOneMax_ = 1.0 / Tmax_;
            log10Pmin_ = std::log10(Pmin_);
            log10Pmax_ = std::log10(Pmax_);
        }

        /**\brief Construct from components.*/
        inline ChebyshevPolynomialsRR(const realArrayArray &ANM)
            : reactionRateTypes(), Anm_(ANM), Tmin_(ChebyshevRRLimit::temperatureLimit::Tmin),
              Tmax_(ChebyshevRRLimit::temperatureLimit::Tmax),
              Pmin_(ChebyshevRRLimit::pressureLimit::Pmin),
              Pmax_(ChebyshevRRLimit::pressureLimit::Pmax), thirdBodyEfficiency_(nullptr) {
            powTminusOneMin_ = 1.0 / Tmin_;
            powTminusOneMax_ = 1.0 / Tmax_;
            log10Pmin_ = std::log10(Pmin_);
            log10Pmax_ = std::log10(Pmax_);
        }

        /**\brief Construct from components.*/
        inline ChebyshevPolynomialsRR(const realArrayArray &ANM, const thirdBodyEfficiency &tbe)
            : reactionRateTypes(), Anm_(ANM), Tmin_(ChebyshevRRLimit::temperatureLimit::Tmin),
              Tmax_(ChebyshevRRLimit::temperatureLimit::Tmax),
              Pmin_(ChebyshevRRLimit::pressureLimit::Pmin),
              Pmax_(ChebyshevRRLimit::pressureLimit::Pmax), thirdBodyEfficiency_(nullptr) {
            powTminusOneMin_ = 1.0 / Tmin_;
            powTminusOneMax_ = 1.0 / Tmax_;
            log10Pmin_ = std::log10(Pmin_);
            log10Pmax_ = std::log10(Pmax_);
            thirdBodyEfficiency_.reset(new thirdBodyEfficiency(tbe));
        }

        /**
         * \brief Construct from components
         * \param[in] defaultTorP - false: using default pressure limits
         * \param[in] min - Tmin
         * \param[in] max - Tmax
         *       Vise verse.
         */
        ChebyshevPolynomialsRR(const realArrayArray &ANM, const real min, const real max,
                               const bool defaultTorP);

        /**
         * \brief Construct from components
         * \param[in] defaultTorP - false: using default pressure limits
         * \param[in] min - Tmin
         * \param[in] max - Tmax
         *       Vise verse.
         */
        ChebyshevPolynomialsRR(const realArrayArray &ANM, const real min, const real max,
                               const bool defaultTorP, const thirdBodyEfficiency &tbe);

        ChebyshevPolynomialsRR(const realArrayArray &ANM, const real Tmin, const real Tmax,
                               const real Pmin, const real Pmax);

        ChebyshevPolynomialsRR(const realArrayArray &ANM, const real Tmin, const real Tmax,
                               const real Pmin, const real Pmax, const thirdBodyEfficiency &tbe);

        ChebyshevPolynomialsRR(const realArray &ANM, const realArray &TMinMax,
                               const realArray &PMinMax, const speciesList &spt,
                               const realArray &efficiences);

        inline ChebyshevPolynomialsRR(const ChebyshevPolynomialsRR &other)
            : reactionRateTypes(other), Anm_(other.Anm_), Tmin_(other.Tmin_),
              powTminusOneMin_(other.powTminusOneMin_), Tmax_(other.Tmax_),
              powTminusOneMax_(other.powTminusOneMax_), Pmin_(other.Pmin_),
              log10Pmin_(other.log10Pmin_), Pmax_(other.Pmax_), log10Pmax_(other.log10Pmax_),
              thirdBodyEfficiency_(nullptr) {
            if (other.thirdBodyEfficiency_) {
                thirdBodyEfficiency_.reset(new thirdBodyEfficiency(*(other.thirdBodyEfficiency_)));
            }
        }

        ChebyshevPolynomialsRR(const speciesList &sp, const controller &cont);

        inline virtual ~ChebyshevPolynomialsRR() noexcept {}

        hur_nodiscard virtual inline uniquePtr<reactionRateTypes> clone() const {
            return uniquePtr<reactionRateTypes>(new ChebyshevPolynomialsRR(*this));
        }

        hur_nodiscard virtual inline std::string type() const noexcept {
            return std::string("usingChebyshevPolynomialsReactionRate");
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

        /**\brief The n x m coefficients.*/
        hur_nodiscard inline realArrayArray &Anm() noexcept { return Anm_; }
        hur_nodiscard inline const realArrayArray &Anm() const noexcept { return Anm_; }
        hur_nodiscard inline realArray &tmpAnm() noexcept { return tmpAnm_; }

        void restructFromTmpAnm();

        hur_nodiscard inline real Tmin() const noexcept { return Tmin_; }
        hur_nodiscard inline real &Tmin() noexcept { return Tmin_; }

        hur_nodiscard inline real Tmax() const noexcept { return Tmax_; }
        hur_nodiscard inline real &Tmax() noexcept { return Tmax_; }

        hur_nodiscard inline real Pmin() const noexcept { return Pmin_; }
        hur_nodiscard inline real &Pmin() noexcept { return Pmin_; }

        hur_nodiscard inline real Pmax() const noexcept { return Pmax_; }
        hur_nodiscard inline real &Pmax() noexcept { return Pmax_; }

        void resetThirdBodyEff(const thirdBodyEfficiency &tbe);

        hur_nodiscard inline virtual bool isModefiedByThirdBody() const noexcept {
            return !thirdBodyEfficiency_.isNull();
        }
        hur_nodiscard inline virtual bool isPressureDenpendent() const noexcept { return true; }
    };
} // namespace OpenHurricane
