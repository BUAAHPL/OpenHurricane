/*!
 * \file combustionModel.hpp
 * \brief Header of combustion model.
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
#include "chemistrySource.hpp"
#include "turbulenceModel.hpp"

namespace OpenHurricane {
    class timeMarching;
    class spatialScheme;

    class combustionModel {
    public:
        enum class chemistryOption : short {
            coupled =
                0, //!<  directly use chemical source terms, i.e. fully coupled. Valid for laminar finite rate, PaSR and EDC as well as ED model.
            strangSplitted =
                1, //!< Strang operator-splitting method. Valid for laminar finite rate model.
            integrated =
                2 //!< Approximate the reaction rate with integrated one. Valid for laminar finite rate, PaSR and EDC model.
        };

    protected:
        flowModel &flows_;

        uniquePtr<chemistrySource> chemistryPtr_;

        turbulenceModel &turbulence_;

        integer maxSubStep_;

        stringList fuelName_;
        stringList oxygenName_;
        stringList productionName_;

        /** \brief The minimum factor to reduce the time step. Default is 0.5. */
        real minDtFactor_;

        mutable integerList *fuelOxygenProduIdPtr_;

        /**
         * \brief chemistry option. Default is coupled.
         */
        chemistryOption chemOptions_;

        /**
         * \brief The time-step for calculate chemical source in integrated scheme.
         */
        mutable realArray *tauIntPtr_;

        /**
         * \brief The time-step for calculate chemical source in integrated scheme.
         */
        hur_nodiscard inline realArray &tauInt() noexcept;

        /**
         * \brief The mixture fraction Z.
         */
        mutable cellRealArray *mixZPtr_;

        realArray fuelMassFraction_;
        realArray oxygenMassFraction_;

        integerList fuelId_;
        integerList oxygenId_;
        real stoichiometricRatio_;

    public:
        static std::map<std::string, chemistryOption> validOptions;

        static void addValidOptions(const std::string &opt, const chemistryOption param);

        //! \cond internalClass
        class initValidTables {
        public:
            initValidTables();
        };
        //! \endcond

        declareClassNames;
        declareObjFty(combustionModel, controller,
                      (flowModel & flows, const controller &cont, turbulenceModel &turb),
                      (flows, cont, turb));

        // Constructors.

        /*!\brief Disallow null constructor.*/
        combustionModel() = delete;

        /*!\brief Construct from flow, controller and turbulence model.*/
        combustionModel(flowModel &flows, const controller &cont, turbulenceModel &turb);

        /*!\brief Disallow copy constructor.*/
        combustionModel(const combustionModel &) = delete;

        static uniquePtr<combustionModel> creator(flowModel &flows, const controller &cont,
                                                  turbulenceModel &turb);

        /*!\brief Destructor.*/
        virtual ~combustionModel() noexcept {
            HurDelete(fuelOxygenProduIdPtr_);
            HurDelete(tauIntPtr_);
            HurDelete(mixZPtr_);
        }

        /*!\brief Return reaction table.*/
        hur_nodiscard inline reactionList &reactions();

        /*!\brief Return reaction table.*/
        hur_nodiscard inline const reactionList &reactions() const;

        /*!\brief Return chemistry source term method.*/
        inline chemistrySource &chemistry() noexcept;

        /*!\brief Return species list.*/
        inline speciesList &species();
        /*!\brief Return species list.*/
        inline const speciesList &species() const;

        /*!\brief Return the mixtures.*/
        inline mixture &mixtures();

        /*!\brief Return the mixtures.*/
        hur_nodiscard inline const mixture &mixtures() const;

        /**
         * \brief The name list of fule.
         */
        hur_nodiscard inline const stringList &fuelName() const noexcept;

        /**
         * \brief The name list of oxygen.
         */
        hur_nodiscard inline const stringList &oxygenName() const noexcept;

        /**
         * \brief The name list of production.
         */
        hur_nodiscard inline const stringList &productionName() const noexcept;

        /**
         * \brief Explicit chemistry source.
         */
        virtual void expChemistrySource(realArray &dt, const bool isModifiedDt = false) = 0;

        /**
         * \brief Implicit chemistry source.
         */
        virtual void impChemistrySource(realArray &dt, const bool isModifiedDt = false) = 0;

        virtual void getChemistrySource() = 0;
        virtual void getImpChemistrySource(realArray &dt, const bool isModifiedDt = false) = 0;

        /**
         * \brief Full point-implicit for NS equations.
         *        The index of species in Jacobian matrix must be continuous.
         * \param[in,out] dt - Timestep
         * \param[out] Jac - The chemical source terms Jacobian based on conservative variables
         * \param[in] rhoId - The index of density rho in Jac
         * \param[in] rhouId - The index of rhou in Jac (rhovId=rhouid+1, rhowId=rhovId+1)
         * \param[in] rhoEId - The index of rhoE in Jac
         * \param[in] rhoYi0Id - The index of rhoYi0 in Jac (rhoYi1Id=rhoYi0Id+1,...)
         */
        virtual void fullPointImpChemistrySource(realArray &dt, cellRealSquareMatrixArray &Jac,
                                                 const integer rhoId, const integer rhouId,
                                                 const integer rhoEId, const integer rhoYi0Id);

        /**
         * \brief Return heat release rate [J/(m^3 s)].
         */
        virtual realArray heatReleaseRate() = 0;

        /*!\brief Solve a constant volume reactor for all internal field.
         * The species mass fraction field and the temperature field will be changed.
         */
        virtual void constVolReactor(const realArray &dt, const real dtFactor,
                                     const bool firstCall = true) {}

        /*!\brief Disallow bitwise assignment.*/
        void operator=(const combustionModel &) = delete;

        /**
         * \brief Damkohler number.
         */
        virtual realArray calcDamkohler();

        /*!\brief Chemical source terms.*/
        virtual realArrayArray omegai() const = 0;

        /**
         * \brief Flame index defined by using the approach of Lacaze.
         * Ref [1] LACAZE G, RICHARDSON E, POINSOT T. Large eddy simulation of spark ignition in a turbulent methane jet[J/OL]. Combustion and Flame, 2009, 156(10): 1993ï¿½C2009.
         *  http://dx.doi.org/10.1016/j.combustflame.2009.05.006. DOI:10.1016/j.combustflame.2009.05.006.
         *
         */
        hur_nodiscard realArray flameIndex() const;

        /*!
         * \brief The characteristic chemical time scale using the approach based on the forward reaction rates.
         *
         */
        inline virtual hur_nodiscard realArray tcFRR();

        /*!
         * \brief The characteristic chemical time scale using the approach based on the forward reaction rates.
         *
         */
        virtual hur_nodiscard void tcFRR(realArray &tc);

        /*!
         * \brief The characteristic chemical time scale using the approach based on the species formation rates.
         *
         */
        inline virtual hur_nodiscard realArray tcSFR();

        /*!
         * \brief The characteristic chemical time scale using the approach based on the species formation rates.
         *
         */
        virtual hur_nodiscard void tcSFR(realArray &tc) const;

        /*!
         * \brief The characteristic chemical time scale using the approach based on the given species production rates.
         *
         */
        virtual hur_nodiscard realArray tcGSPR(const integer isp);

        /*!
         * \brief The characteristic chemical time scale using the approach based on the fuel, oxygen and species production.
         *
         */
        inline virtual hur_nodiscard realArray tcGSPR();

        /*!
         * \brief The characteristic chemical time scale using the approach based on the given species.
         * \param[in] maxOrMin - True = maximum; False = minimum
         */
        virtual void tcGSPR(realArray &tc, const integerList &spid,
                            const bool maxOrMin = true) const;

        /*!
         * \brief The characteristic chemical time scale using the approach based on the fuel, oxygen and species production.
         *
         */
        virtual hur_nodiscard void tcGSPR(realArray &tc) const;

        /*!
         * \brief The characteristic chemical time scale using the approach based on the Jacobian diagonal terms.
         *
         */
        inline virtual hur_nodiscard realArray tcJacDT();

        /*!
         * \brief The characteristic chemical time scale using the approach based on the Jacobian diagonal terms.
         *
         */
        virtual hur_nodiscard void tcJacDT(realArray &tc);

        hur_nodiscard integerList fuelOxygenProduId() const;

        /**
         * \brief Calculate the chemical source terms of species transport equations.
         *        The index of species in Jacobian matrix must be continuous.
         * \param[in,out] times - Time marching method
         * \param[in] rhoId - The index of density rho in Jac
         * \param[in] rhouId - The index of rhou in Jac (rhovId=rhouid+1, rhowId=rhovId+1)
         * \param[in] rhoEId - The index of rhoE in Jac
         * \param[in] rhoYi0Id - The index of rhoYi0 in Jac (rhoYi1Id=rhoYi0Id+1,...)
         */
        virtual void evaluateSource(timeMarching &times, const integer rhoId, const integer rhouId,
                                    const integer rhoEId, const integer rhoYi0Id);

        virtual void integratedChemistrySource();

        /**
         * \brief Calculate the chemical source terms in strang operator-splitting method.
         * \param[in,out] times - Time marching method
         */
        virtual void strangSplittedChemistrySource(timeMarching &times);

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept;

        hur_nodiscard inline bool isStrangSplitted() const noexcept;
        hur_nodiscard inline bool isCoupled() const noexcept;
        hur_nodiscard inline bool isIntegrated() const noexcept;

    protected:
        /**
         * \brief Calculate the chemical source terms of species transport equations in coupled scheme.
         *        The index of species in Jacobian matrix must be continuous.
         * \param[in,out] times - Time marching method
         * \param[in] rhoId - The index of density rho in Jac
         * \param[in] rhouId - The index of rhou in Jac (rhovId=rhouid+1, rhowId=rhovId+1)
         * \param[in] rhoEId - The index of rhoE in Jac
         * \param[in] rhoYi0Id - The index of rhoYi0 in Jac (rhoYi1Id=rhoYi0Id+1,...)
         */
        virtual void chemistrySourceCoupled(timeMarching &times, const integer rhoId,
                                            const integer rhouId, const integer rhoEId,
                                            const integer rhoYi0Id);

        /**
         * \brief Caculate chemical source terms explicitly in coupled scheme.
         */
        virtual void expChemistrySourceCoupled() = 0;

        /**
         * \brief Caculate chemical source terms implicitly with diagonal Jacobian matrices in coupled scheme.
         */
        virtual void diagImpChemistrySourceCoupled() = 0;

        virtual void fullImpChemistrySourceCoupled(cellRealSquareMatrixArray &Jac,
                                                   const integer rhoId, const integer rhouId,
                                                   const integer rhoEId,
                                                   const integer rhoYi0Id) = 0;

        /**
         * \brief Caculate chemical source terms in integrated scheme.
         */
        virtual void chemistrySourceIntegrated() = 0;

        /**
         * \brief Get chemical source terms in integrated scheme to species equations.
         */
        virtual void getIntegratedChemSource();

        /**
         * \brief Caculate chemical source terms in strangSplitted scheme.
         */
        virtual void chemistrySourceStrangSplitted(const realArray &dt, const real dtFactor) = 0;

    public: // About mixture
        bool calcMixtureFraction();

        hur_nodiscard inline cellRealArray &mixZ();
        hur_nodiscard inline const cellRealArray &mixZ() const;

        realArray mixtureDisspationRate(const spatialScheme &sps);
    };
} // namespace OpenHurricane

#include "combustionModel.inl"