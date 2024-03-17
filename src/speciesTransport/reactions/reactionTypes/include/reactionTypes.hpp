/*!
 * \file reactionTypes.hpp
 * \brief Header of base class of reaction.
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
#include "objectFactory.hpp"
#include "reactionRateTypes.hpp"

#ifdef CUDA_PARALLEL
#include "cuChemInterface.hpp"
#endif // CUDA_PARALLEL

namespace OpenHurricane {
    // Forward declaration of reaction table class.
    class reactionList;

    /*!\brief The base class of reaction.*/
    class reactionTypes {

    public:
        enum class types : short {
            irreversible = 0,
            reversible = 1,
            nonEqReversible = 2,
            unknown = 3
        };

        struct reactionSpeciesCoeffs {

            /*!\brief Index of species.*/
            integer index_;

            /*!\brief stoichiometry coefficients.*/
            real stoichCoeff_;

            /*!\brief Reaction order.*/
            real order_;
        };

    protected:
        // Private data

        /*!\brief Name of the reaction.*/
        const string name_;

        /*!\brief Const reference to the species list.*/
        const speciesList &species_;

        /*!\brief Const reference to the Gibbs free energy of species at standard pressure.*/
        const realArray &G0_;

        /*!\brief Const reference to the Gibbs free energy of species at standard pressure.*/
        const realArray &DG0DT_;

        /*!\brief Index of reaction.*/
        integer index_;

        /*!\brief Forward species stoichiometry coeffs of the reaction.*/
        List<reactionSpeciesCoeffs> forwardCoeffs_;

        /*!\brief Backward species stoichiometry coeffs of the reaction.*/
        List<reactionSpeciesCoeffs> backwardCoeffs_;

        /*!\brief backwardCoeffs_ - forwardCoeffs_.*/
        List<reactionSpeciesCoeffs> bfCoeffs_;

        /*!\brief The summation of stoichiometry coeffs.*/
        real sumStoi_;

    protected:
        types reacType_;

    public:
        declareClassName(reactionTypes);
        declareObjFty(reactionTypes, controller,
                      (const speciesList &sp, const reactionList &rt,
                       const controller &reactionControl),
                      (sp, rt, reactionControl));

        // Constructors

        /*!\brief Construct from components.*/
        reactionTypes(const string &names, const speciesList &specie, const reactionList &rt,
                      const integer id, const List<reactionSpeciesCoeffs> &forwardCoe,
                      const List<reactionSpeciesCoeffs> &backwardCoe);

        /*!\brief Construct from components.*/
        reactionTypes(const string &names, const speciesList &specie, const reactionList &rt,
                      const integer id);

        /*!\brief Construct from controller*/
        reactionTypes(const speciesList &sp, const reactionList &rt, const controller &cont);

        /*!\brief Construct as copy.*/
        inline reactionTypes(const reactionTypes &rs)
            : name_(rs.name_), species_(rs.species_), G0_(rs.G0_), DG0DT_(rs.DG0DT_),
              index_(rs.index_), forwardCoeffs_(rs.forwardCoeffs_),
              backwardCoeffs_(rs.backwardCoeffs_), bfCoeffs_(rs.bfCoeffs_), sumStoi_(rs.sumStoi_),
              reacType_(rs.reacType_) {}

        /*!\brief Return a clone.*/
        hur_nodiscard virtual inline uniquePtr<reactionTypes> clone() const {
            return uniquePtr<reactionTypes>(new reactionTypes(*this));
        }

        hur_nodiscard static uniquePtr<reactionTypes>
        creator(const speciesList &sp, const reactionList &rt, const controller &cont);

        /*!\brief Destructor.*/
        virtual inline ~reactionTypes() noexcept {}

        // Member Functions

        /*!\brief Return const access to the reaction name.*/
        hur_nodiscard inline const string &name() const noexcept { return name_; }

        /*!\brief Return const reference to species list.*/
        hur_nodiscard inline const speciesList &species() const noexcept { return species_; }

        /*!\brief Forward species stoichiometry coeffs of the reaction.*/
        hur_nodiscard inline const List<reactionSpeciesCoeffs> &forwardCoeffs() const noexcept {
            return forwardCoeffs_;
        }

        /*!\brief Forward species stoichiometry coeffs of the reaction.*/
        hur_nodiscard inline List<reactionSpeciesCoeffs> &forwardCoeffs() noexcept {
            return forwardCoeffs_;
        }

        /*!\brief Backward species stoichiometry coeffs of the reaction.*/
        hur_nodiscard inline const List<reactionSpeciesCoeffs> &backwardCoeffs() const noexcept {
            return backwardCoeffs_;
        }

        /*!\brief Backward species stoichiometry coeffs of the reaction.*/
        hur_nodiscard inline List<reactionSpeciesCoeffs> &backwardCoeffs() noexcept {
            return backwardCoeffs_;
        }

        hur_nodiscard inline const List<reactionSpeciesCoeffs> &bfCoeffs() const noexcept {
            return bfCoeffs_;
        }

        hur_nodiscard inline List<reactionSpeciesCoeffs> &bfCoeffs() noexcept { return bfCoeffs_; }
        void setBFCoeffs();

        /*!\brief Index of reaction.*/
        hur_nodiscard inline integer index() const noexcept { return index_; }

        // Reaction rate constant

        /*!\brief Forward rate constant.*/
        hur_nodiscard virtual inline real kf(const real p, const real T, const realArray &c) const {
            return real(0.0);
        }

        /*!\brief Reverse rate constant from the given forward constant.*/
        hur_nodiscard virtual inline real kr(const real kfwd, const real p, const real T,
                                             const realArray &c) const {
            return real(0.0);
        }

        /*!\brief Reverse rate constant.*/
        hur_nodiscard virtual inline real kr(const real p, const real T, const realArray &c) const {
            return real(0.0);
        }

        /*!\brief Equilibrium constant i.t.o. partial pressures*/
        hur_nodiscard real Kp(const real T) const;

        /*!\brief Temperature derivative of equilibrium constant [] i.t.o. partial pressures*/
        hur_nodiscard real DKpDT(const real T) const;

        /*!\brief Equilibrium constant i.t.o. molar concentration.*/
        hur_nodiscard real Kc(const real T) const;

        /*!\brief Temperature derivative of equilibrium constant i.t.o. molar concentration.*/
        hur_nodiscard real DKcDT(const real T) const;

        /*!\brief Net reaction rate.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \param[out] qf - forward reaction rate [kmol/m^3/s]
         * \param[out] qr - reverse reaction rate [kmol/m^3/s]
         * \return Return the net reaction rate i.e. (qf - qr)
         */
        real q(const real p, const real T, const realArray &c, real &qf, real &qr) const;

        /*!\brief Net reaction rate.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \param[out] qf - forward reaction rate [kmol/m^3/s]
         * \param[out] qr - reverse reaction rate [kmol/m^3/s]
         * \return Return the net reaction rate i.e. (qf - qr)
         */
        real q(const real p, const real T, const realArray &c, real &qf, real &qr, real &kf,
               real &kr) const;

        /*!\brief Return Net reaction rate and add to each species source terms.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \param[in] dcidt - species source terms [kmol/m^3/s]
         * \return Return the net reaction rate.
         */
        real q(const real p, const real T, const realArray &c, realArray &dcidt) const;

        /*!\brief Return Net reaction rate and add to each species source terms.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \param[in] dcidt - species source terms [kmol/m^3/s]
         * \return Return the net reaction rate.
         */
        real q(const real p, const real T, const realArray &c, realArray &dcidt, real &qf, real &qr,
               real &kf, real &kr) const;

        /*!\brief Return Net reaction rate and add to each species source terms.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \param[in] dcidt - species source terms [kmol/m^3/s]
         * \return Return the net reaction rate.
         */
        real q(const real p, const real T, const realArray &c, const bool reduced,
               const integerArray &fullToSimplified, realArray &dcidt) const;

        /*!\brief Return Net reaction rate and add to each species source terms.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \param[in] dcidt - species source terms [kmol/m^3/s]
         * \return Return the net reaction rate.
         */
        real q(const real p, const real T, const realArray &c, const bool reduced,
               const integerArray &fullToSimplified, realArray &dcidt, real &qf, real &qr, real &kf,
               real &kr) const;

        /*!\brief Return the temperature derivative of forward rate.
         * \param[in] kf - forward rate
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \return Return the temperature derivative of forward rate.
         */
        hur_nodiscard virtual real DkfDT(const real kf, const real p, const real T,
                                         const realArray &c) const {
            checkWarning("This function must be overwritten in the derived class.");
            return real(0);
        }

        /*!\brief Return the temperature derivative of forward rate.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \return Return the temperature derivative of forward rate.
         */
        hur_nodiscard virtual real DkfDT(const real p, const real T, const realArray &c) const {
            checkWarning("This function must be overwritten in the derived class.");
            return real(0);
        }

        /*!\brief Return the temperature derivative of backward rate.
         * \param[in] kr - backward rate
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \return Return the temperature derivative of backward rate.
         */
        hur_nodiscard virtual real DkrDT(const real kr, const real p, const real T,
                                         const realArray &c, const real dkfdT) const {
            checkWarning("This function must be overwritten in the derived class.");
            return real(0);
        }

        /** \brief The partial derivatives of third-body terms with respect to the temperature.*/
        hur_nodiscard virtual real DGGDT(const real p, const real T, const realArray &c) const {
            checkWarning("This function must be overwritten in the derived class.");
            return real(0);
        }

        /** \brief The partial derivatives of third-body terms with respect to the species molar concentrations.*/
        hur_nodiscard virtual void DGGDc(const real p, const real T, const realArray &c,
                                         realArray &dgdci) const {
            checkWarning("This function must be overwritten in the derived class.");
        }

        /*!\brief Return the pressure derivative of forward rate.
         * \param[in] kf - forward rate
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \return Return the temperature derivative of forward rate.
         */
        hur_nodiscard virtual real DkfDp(const real kf, const real p, const real T,
                                         const realArray &c) const {
            checkWarning("This function must be overwritten in the derived class.");
            return real(0);
        }

        /*!\brief Return the pressure derivative of forward rate.
         * \param[in] kr - forward rate
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \return Return the temperature derivative of forward rate.
         */
        hur_nodiscard virtual real DkrDp(const real kr, const real p, const real T,
                                         const realArray &c, const real dkfdp) const {
            checkWarning("This function must be overwritten in the derived class.");
            return real(0);
        }

        void DqDci(const real kf, const real kr, const real qfr, const real p, const real T,
                   const realArray &c, realSquareMatrix &Jac) const;

        /*!\brief Calculate the diagonal of Jacobian of chemical source terms d(Ri)/d(rhoi).
         * \param[in] kf - forward reaction rate.
         * \param[in] kr - reverse reaction rate.
         * \param[in] qfr - net reaction rate.
         * \param[in] p - pressure
         * \param[in] T -temperature
         * \param[in] c - molar concentrations (Should be full size of all species).
         * \param[out] diagJac - The diagonal terms of chemical source Jacobian.
         * \param[in] calcLastSpc - Should treat the last species as independent parameter.
         */
        void DqDci(const real kf, const real kr, const real qfr, const real p, const real T,
                   const realArray &c, realArray &diagJac, const bool calcLastSpc) const;

        /*!\brief Return Net reaction rate and add to each species source terms.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \param[in] dcidt - species source terms [kmol/m^3/s]
         * \return Return the net reaction rate.
         */
        real q(const real p, const real T, const realArray &c, realArray &dcidt, realArray &diagJac,
               const bool calcLastSpc, real &kf, real &kr) const;

        void DqDci(const real kf, const real kr, const real qfr, const real p, const real T,
                   const realArray &c, const bool reduced, const integerArray &fullToSimplified,
                   realSquareMatrix &Jac) const;

        /*!\brief Calculate the diagonal of Jacobian of chemical source terms d(Ri)/d(rhoi).
         * \param[in] kf - forward reaction rate.
         * \param[in] kr - reverse reaction rate.
         * \param[in] qfr - net reaction rate.
         * \param[in] p - pressure
         * \param[in] T -temperature
         * \param[in] c - molar concentrations (Should be full size of all species).
         * \param[in] reduced - If the mechanism reduced?
         * \param[in] fullToSimplified - the map from full species list to simplified species list
         * \param[out] diagJac - The diagonal terms of chemical source Jacobian.
         * \param[in] calcLastSpc - Should treat the last species as independent parameter.
         */
        void DqDci(const real kf, const real kr, const real qfr, const real p, const real T,
                   const realArray &c, const bool reduced, const integerArray &fullToSimplified,
                   realArray &diagJac, const bool calcLastSpc) const;

        /*!\brief Return Net reaction rate and add to each species source terms.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration of species [kmol/m^3]
         * \param[in] dcidt - species source terms [kmol/m^3/s]
         * \return Return the net reaction rate.
         */
        real q(const real p, const real T, const realArray &c, const bool reduced,
               const integerArray &fullToSimplified, realArray &dcidt, realArray &diagJac,
               const bool calcLastSpc, real &kf, real &kr) const;

        void DqDT(const real kf, const real kr, const real qfr, const real p, const real T,
                  const realArray &c, realSquareMatrix &Jac) const;

        void DqDT(const real kf, const real kr, const real qfr, const real p, const real T,
                  const realArray &c, realArray &dwdT) const;

        void DqDT(const real kf, const real kr, const real qfr, const real p, const real T,
                  const realArray &c, const bool reduced, const integerArray &fullToSimplified,
                  const integer Tid, realSquareMatrix &Jac) const;

        void DqDT(const real kf, const real kr, const real qfr, const real p, const real T,
                  const realArray &c, const bool reduced, const integerArray &fullToSimplified,
                  realArray &dwdT) const;

        void DqDp(const real kf, const real kr, const real qfr, const real p, const real T,
                  const realArray &c, realSquareMatrix &Jac) const;

        void DqDp(const real kf, const real kr, const real qfr, const real p, const real T,
                  const realArray &c, realArray &dwdp) const;

        void DqDp(const real kf, const real kr, const real qfr, const real p, const real T,
                  const realArray &c, const bool reduced, const integerArray &fullToSimplified,
                  const integer pid, realSquareMatrix &Jac) const;

        void DqDp(const real kf, const real kr, const real qfr, const real p, const real T,
                  const realArray &c, const bool reduced, const integerArray &fullToSimplified,
                  realArray &dwdp) const;

        /*!\brief Return true if the reaction contains the third-body terms.*/
        hur_nodiscard virtual bool hasThirdBodyTerms() const noexcept {
            checkWarning("This function must be overwritten in the derived class.");
            return false;
        }

        hur_nodiscard virtual bool isPressureDenpendent() const noexcept {
            checkWarning("This function must be overwritten in the derived class.");
            return false;
        }

        hur_nodiscard inline types reacType() const noexcept { return reacType_; }

        hur_nodiscard inline virtual bool isModefiedByThirdBody() const noexcept {
            LFatal("This function must be overwritten in the derived class");
            return false;
        }

        /**
         * \brief The partial derivatives of third-body terms with respect to the molar concentration of species.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        virtual inline void gamDGamDCi(const real p, const real T, const realList &c,
                                       realList &gdgdci) const {
            LFatal("This function must be overwritten in the derived class");
        }

        /**
         * \brief The partial derivatives of third-body terms with respect to the temperature.
         * \param[in] p - The static pressure [Pa].
         * \param[in] T - The static temperature [K].
         * \param[in] c - The molar concentrations of species [kmol/m^3].
         */
        hur_nodiscard virtual inline real gamDGamDT(const real p, const real T,
                                                    const realList &c) const {
            LFatal("This function must be overwritten in the derived class");
            return 0;
        }

        
#ifdef CUDA_PARALLEL

    protected:
        void cuSetReactionCoeff(const integer reacId, const integer nrc,
                                cuChem::cuChemInterface::reactionCoeffsInterface &reacInt,
                                const List<reactionSpeciesCoeffs> &coef) const;

    public:
        virtual void cuSetReactionType(const integer reacId, const integer nrc,
                                       cuChem::cuChemInterface &reacInt) const;

#endif // CUDA_PARALLEL
    };

} // namespace OpenHurricane