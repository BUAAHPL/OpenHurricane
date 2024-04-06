/*!
 * \file fallOffFunctions.hpp
 * \brief Header of fall-off functions.
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
#include "controller.hpp"
#include "objectFactory.hpp"
#include "real.hpp"
#include "smartPointer.hpp"

namespace OpenHurricane {
    namespace fallOffFunctions {

        /**
         * \brief The base class of fall-off functions.
         */
        class fallOffFunction {
        private:
        public:
            declareClassName(fallOffFunction);
            declareObjFty(fallOffFunction, controller, (const controller &cont), (cont));

            fallOffFunction() = default;

            inline fallOffFunction(const fallOffFunction &fof) {}

            inline fallOffFunction(const controller &cont) {}

            hur_nodiscard static uniquePtr<fallOffFunction> creator(const controller &cont);

            virtual inline ~fallOffFunction() noexcept {}

            hur_nodiscard virtual inline uniquePtr<fallOffFunction> clone() const {
                return uniquePtr<fallOffFunction>(nullptr);
            }

            virtual hur_nodiscard inline std::string type() const noexcept = 0;

            /**
             * \brief The operator for calculating blending function.
             * \param[in] T - The static temperature.
             * \param[in] Pr - The reduced pressure.
             * \return The blending function F.
             */
            hur_nodiscard virtual inline real F(const real T, const real Pr) const = 0;

            hur_nodiscard virtual inline real DFDci(const real T, const real Pr, const real F,
                                                    const real DPrDci) const = 0;

            hur_nodiscard virtual inline real DFDT(const real T, const real Pr, const real F,
                                                   const real DPrDT) const = 0;

            hur_nodiscard inline virtual bool isLindemann() const noexcept { return false; }
            hur_nodiscard inline virtual bool isSRI() const noexcept { return false; }
            hur_nodiscard inline virtual bool isTroe() const noexcept { return false; }
        };

        /**
         * \brief The class of Lindemann fall-off function.
         */
        class Lindemann : public fallOffFunction {
        public:
            declareClassName(Lindemann);
            inline Lindemann() : fallOffFunction() {}

            inline Lindemann(const Lindemann &other) : fallOffFunction(other) {}

            inline Lindemann(const controller &cont) : fallOffFunction(cont) {}

            virtual inline ~Lindemann() noexcept {}

            hur_nodiscard virtual inline uniquePtr<fallOffFunction> clone() const {
                return uniquePtr<fallOffFunction>(new Lindemann(*this));
            }

            virtual hur_nodiscard inline std::string type() const noexcept {
                return std::string("Lindemann");
            }

            /**
             * \brief The operator for calculating blending function.
             * \param[in] T - The static temperature.
             * \param[in] Pr - The reduced pressure.
             * \return The blending function F.
             */
            hur_nodiscard virtual inline real F(const real T, const real Pr) const {
                return real(1.0);
            }

            hur_nodiscard virtual inline real DFDci(const real T, const real Pr, const real F,
                                                    const real DPrDci) const {
                return real(0);
            }

            hur_nodiscard virtual inline real DFDT(const real T, const real Pr, const real F,
                                                   const real DPrDT) const {
                return real(0);
            }

            hur_nodiscard inline virtual bool isLindemann() const noexcept { return true; }
        };

        /**
         * \brief The class of SRI fall-off function.
         */
        class SRI : public fallOffFunction {
        private:
            // Private data

            real a_;
            real b_;
            real c_;
            real d_;
            real e_;

        public:
            declareClassName(SRI);

            /**\brief Construct null.*/
            inline SRI() : fallOffFunction(), a_(0.0), b_(0.0), c_(0.0), d_(1.0), e_(0.0) {}

            /**\brief Construct from components with three parameters.*/
            inline SRI(const real a, const real b, const real c)
                : fallOffFunction(), a_(a), b_(b), c_(c), d_(1.0), e_(0.0) {}

            /**\brief Construct from components with five parameters.*/
            inline SRI(const real a, const real b, const real c, const real d, const real e)
                : fallOffFunction(), a_(a), b_(b), c_(c), d_(d), e_(e) {}

            inline SRI(const SRI &other)
                : fallOffFunction(other), a_(other.a_), b_(other.b_), c_(other.c_), d_(other.d_),
                  e_(other.e_) {}

            /*!\brief Construct from controller.*/
            inline SRI(const controller &cont)
                : fallOffFunction(cont), a_(cont.findOrDefault<real>("a", real(0.0))),
                  b_(cont.findOrDefault<real>("b", real(0.0))),
                  c_(cont.findOrDefault<real>("e", real(0.0))),
                  d_(cont.findOrDefault<real>("d", real(1.0))),
                  e_(cont.findOrDefault<real>("e", real(0.0))) {}

            virtual inline ~SRI() noexcept {}

            hur_nodiscard virtual inline uniquePtr<fallOffFunction> clone() const {
                return uniquePtr<fallOffFunction>(new SRI(*this));
            }

            virtual hur_nodiscard inline std::string type() const noexcept {
                return std::string("SRI");
            }

            /**
             * \brief The operator for calculating blending function.
             * \param[in] T - The static temperature.
             * \param[in] Pr - The reduced pressure.
             * \return The blending function F.
             */
            hur_nodiscard virtual real F(const real T, const real Pr) const;

            hur_nodiscard virtual real DFDci(const real T, const real Pr, const real F,
                                             const real DPrDci) const;

            hur_nodiscard virtual real DFDT(const real T, const real Pr, const real F,
                                            const real DPrDT) const;

            hur_nodiscard inline real a() const noexcept { return a_; }
            hur_nodiscard inline real b() const noexcept { return b_; }
            hur_nodiscard inline real c() const noexcept { return c_; }
            hur_nodiscard inline real d() const noexcept { return d_; }
            hur_nodiscard inline real e() const noexcept { return e_; }
            hur_nodiscard inline virtual bool isSRI() const noexcept { return true; }
        };

        /**
         * \brief The class of Troe fall-off function.
         */
        class Troe : public fallOffFunction {
        private:
            // Private data

            real alpha_;
            real Tsss_, Ts_, Tss_;

            bool isTssOmitted_;

        public:
            declareClassName(Troe);
            inline Troe()
                : fallOffFunction(), alpha_(0.0), Tsss_(0.0), Ts_(0.0), Tss_(0.0),
                  isTssOmitted_(false) {}

            /**
             *\brief- Construct from components with three parameters
             *		¦Á, T***, T*.
             */
            inline Troe(const real alpha, const real Tsss, const real Ts)
                : fallOffFunction(), alpha_(alpha), Tsss_(Tsss), Ts_(Ts), Tss_(0.0),
                  isTssOmitted_(true) {}

            /**
             *\brief Construct from components with four parameters
             *		¦Á, T***, T*, T**.
             */
            inline Troe(const real alpha, const real Tsss, const real Ts, const real Tss)
                : fallOffFunction(), alpha_(alpha), Tsss_(Tsss), Ts_(Ts), Tss_(Tss),
                  isTssOmitted_(false) {}

            /*!\brief Construct from controller.*/
            inline Troe(const Troe &other)
                : fallOffFunction(other), alpha_(other.alpha_), Tsss_(other.Tsss_), Ts_(other.Ts_),
                  Tss_(other.Tss_), isTssOmitted_(other.isTssOmitted_) {}

            inline Troe(const controller &cont)
                : fallOffFunction(cont), alpha_(cont.findOrDefault<real>("alpha", real(0.0))),
                  Tsss_(cont.findOrDefault<real>("Tsss", real(0.0))),
                  Ts_(cont.findOrDefault<real>("Ts", real(0.0))), Tss_(0.0), isTssOmitted_(true) {
                if (cont.found("Tss")) {
                    isTssOmitted_ = false;
                    Tss_ = cont.findOrDefault<real>("Tss", real(0.0));
                }
            }

            virtual inline ~Troe() noexcept {}

            hur_nodiscard virtual inline uniquePtr<fallOffFunction> clone() const {
                return uniquePtr<fallOffFunction>(new Troe(*this));
            }

            virtual hur_nodiscard inline std::string type() const noexcept {
                return std::string("Troe");
            }

            /**
             * \brief The operator for calculating blending function.
             * \param[in] T - The static temperature.
             * \param[in] Pr - The reduced pressure.
             * \return The blending function F.
             */
            hur_nodiscard virtual real F(const real T, const real Pr) const;

            hur_nodiscard virtual real DFDci(const real T, const real Pr, const real F,
                                             const real DPrDci) const;

            hur_nodiscard virtual real DFDT(const real T, const real Pr, const real F,
                                            const real DPrDT) const;

            hur_nodiscard inline real alpha() const noexcept { return alpha_; }
            hur_nodiscard inline real Tsss() const noexcept { return Tsss_; }
            hur_nodiscard inline real Ts() const noexcept { return Ts_; }
            hur_nodiscard inline real Tss() const noexcept { return Tss_; }

            hur_nodiscard inline bool isTssOmitted() const noexcept { return isTssOmitted_; }

            hur_nodiscard inline virtual bool isTroe() const noexcept { return true; }
        };
    } // namespace fallOffFunctions
} // namespace OpenHurricane