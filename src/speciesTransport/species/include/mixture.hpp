/*!
 * \file mixture.hpp
 * \brief Header of class of mixture.
 *       The subroutines and functions are in the <i>mixture.cpp</i> file.
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

#include "reactionList.hpp"
#include "thermoList.hpp"
#include "transportList.hpp"

#ifdef CUDA_PARALLEL
#include "cuChemInterface.hpp"
#include "transportListCUDA.hpp"
#endif // CUDA_PARALLEL

namespace OpenHurricane {

    /*!\brief The basic class of  mixture.*/
    class basicMixture {
    protected:
        /*!\brief Is a singular species?*/
        bool isSingular_;

        const runtimeMesh &mesh_;

        /*!\brief The table of species.*/
        speciesList species_;

        /*!\brief The mass fractions of species.*/
        PtrList<cellRealArray> Yi_;

        /*!\brief The mass diffusion coefficient of species.*/
        PtrList<cellRealArray> Dim_;

        /*!\brief The absolute enthalpy of species.*/
        PtrList<cellRealArray> hi_;

    public:
        /*!\brief Construct from mesh and species size*/
        basicMixture(const runtimeMesh &mesh, const integer speciesSize);

        /*!\brief Construct from mesh and species size*/
        basicMixture(const runtimeMesh &mesh, const speciesList &species);

        /*!\brief Destructor.*/
        virtual ~basicMixture() noexcept {}

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept { return mesh_; }

        /*!\brief Is a singular species?*/
        hur_nodiscard inline bool isSingular() const noexcept { return isSingular_; }

        /*!\brief Return const-reference to the table of species.*/
        hur_nodiscard inline const speciesList &species() const noexcept { return species_; }

        /*!\brief Return access to the table of species.*/
        hur_nodiscard inline speciesList &species() noexcept { return species_; }

        /*!\brief Return const-reference to the mass fractions of species.*/
        hur_nodiscard inline const PtrList<cellRealArray> &Yi() const noexcept { return Yi_; }

        /*!\brief Return access to the mass fractions of species.*/
        hur_nodiscard inline PtrList<cellRealArray> &Yi() noexcept { return Yi_; }

        /*!\brief Return const-reference to the mass fractions of species i.*/
        hur_nodiscard inline const cellRealArray &Yi(const integer spId) const noexcept {
            return Yi_[spId];
        }

        /*!\brief Return access to the mass fractions of species i.*/
        hur_nodiscard inline cellRealArray &Yi(const integer spId) noexcept { return Yi_[spId]; }

        /*!\brief Return const-reference to the mass fractions of species given by name.*/
        hur_nodiscard inline const cellRealArray &Yi(const string &speciesName) const {
            return Yi_[species_.index(speciesName)];
        }

        /*!\brief Return access to the mass fractions of species given by name.*/
        hur_nodiscard inline cellRealArray &Yi(const string &speciesName) {
            return Yi_[species_.index(speciesName)];
        }

        /*!
         * \brief Return access to the mass fractions of all species of cell id=cellId.
         * Can not change Yi in cellId.
         */
        hur_nodiscard realArray Yic(const integer cellId) const;

        /*!
         * \brief Return access to the mass fractions of all species of cell id=cellId.
         * Can not change Yi in cellId.
         */
        hur_nodiscard realArrayArray Yif(const integer faceZoneId) const;

        /*!\brief Return const-reference to the mass diffusion coefficient of species.*/
        hur_nodiscard inline const PtrList<cellRealArray> &Dim() const noexcept { return Dim_; }

        /*!\brief Return access to the mass diffusion coefficient of species.*/
        hur_nodiscard inline PtrList<cellRealArray> &Dim() noexcept { return Dim_; }

        /*!\brief Return const-reference to the mass diffusion coefficient of species i.*/
        hur_nodiscard inline const cellRealArray &Dim(const integer i) const noexcept {
            return Dim_[i];
        }

        /*!\brief Return access to the mass diffusion coefficient of species i.*/
        hur_nodiscard inline cellRealArray &Dim(const integer i) noexcept { return Dim_[i]; }

        /*!\brief Return const-reference to the absolute enthalpy of species.*/
        hur_nodiscard inline const PtrList<cellRealArray> &hi() const noexcept { return hi_; }

        /*!\brief Return access to the absolute enthalpy of species.*/
        hur_nodiscard inline PtrList<cellRealArray> &hi() noexcept { return hi_; }

        /*!\brief Return const-reference to the absolute enthalpy of species i.*/
        hur_nodiscard inline const cellRealArray &hi(const integer i) const noexcept {
            return hi_[i];
        }

        /*!\brief Return access to the absolute enthalpy of species i.*/
        hur_nodiscard inline cellRealArray &hi(const integer i) noexcept { return hi_[i]; }

        /*!\brief Return the molecular weight of species i.*/
        hur_nodiscard inline real Wi(const integer i) const noexcept { return species_.W(i); }

        /*!\brief Return the molecular weight for every species.*/
        hur_nodiscard realArray Wi() const;

        /*!
         * \brief Must be called after creating the speciesTable.
         * It is no need for singular specie to call this subroutine.
         */
        void setYi(const runtimeMesh &mesh);

        /*!
         * \brief Must be called after creating the speciesTable.
         * It is no need for singular specie to call this subroutine.
         */
        void setYi(const runtimeMesh &mesh, const realArray &Yi0);

        /**\brief To get the mass fraction of last species and to make sure the sum of mass fraction of all species in each cells is unity.*/
        void lastSpeAndNormalized();

        /*!
         * \brief The Takeno Flame Index.
         * \param[in] fuelName - The name list of fuel
         * \param[in] oxygenName - The name list of oxygen
         * \return The Takeno Flame Index
         * \retval Return a real array
         */
        hur_nodiscard realArray GFO(const stringList &fuelName, const stringList &oxygenName) const;

        /*!
         * \brief The normalizing, Takeno Flame Index.
         * \param[in] fuelName - The name list of fuel
         * \param[in] oxygenName - The name list of oxygen
         * \return The Takeno Flame Index: 1.0 - indicates premixed combustion; -1.0 indicates non-premixed combustion; 0 - indicates no flame
         * \retval Return a real array
         */
        hur_nodiscard realArray nGFO(const stringList &fuelName,
                                     const stringList &oxygenName) const;
    };

    /*!\brief The class of mixture.*/
    class mixture : public basicMixture {
    public:
        enum class mixtureTabulations : short {
            none = 0,
            massDiffusionTable = 1,
            gasPropertiesTable = 2
        };

    private:
        /**
         * \brief The pointer of thermo table.
         */
        uniquePtr<thermoList> thermoPtr_;

        /**
         * \brief The pointer of transport table.
         */
        uniquePtr<transportList> transportPtr_;

        /**
         * \brief The pointer of reaction table.
         */
        uniquePtr<reactionList> reactionPtr_;

#ifdef CUDA_PARALLEL

        /**
         * \brief The pointer of reaction table in CUDA platform.
         */
        uniquePtr<cuChem::cuChemInterface> cuChemInterfacePtr_;
        void cuAllocReaction();
        void cuParsingReaction();
        void cuCorrectReaction();
        void cuAddThermoToCUDA();
        /**
         * \brief The pointer of transport table in CUDA platform.
         */
        mutable uniquePtr<CUDATransport::transportListCUDA> transportCUDAPtr_;

        void makeTransportTableCUDA() const;
        void makeTransportTableCUDA(const cudaStreams &streams) const;

        mutable uniquePtr<cuChem::speciesTable> speciesCUDAPtr_;

        void makeSpeciesCUDA() const;
        void makeSpeciesCUDA(const cudaStreams &streams) const;

#endif // CUDA_PARALLEL

        bool noReaction_;

        mixtureTabulations mixtureTabulation_;

    public:
        // Constructors

        mixture(const runtimeMesh &mesh, const integer speciesSize, const bool noReaction = false);

        mixture(const runtimeMesh &mesh, const speciesList &species, const bool noReaction = false);

        mixture(const runtimeMesh &mesh, const speciesList &species, const controller &cont,
                const bool inviscous = false);

        mixture(const runtimeMesh &mesh, const controller &cont, const bool inviscous = false);

        /*!\brief Destructor.*/
        virtual ~mixture() noexcept;

        // Member functions

        /*!\brief The thermo-table.*/
        hur_nodiscard inline thermoList &thermalTable() {
            if (!thermoPtr_) {
                checkWarning("Attempt to access to null thermo-table pointer");
            }
            return *thermoPtr_;
        }

        /*!\brief The thermo-table.*/
        hur_nodiscard inline const thermoList &thermalTable() const {
            if (!thermoPtr_) {
                checkWarning("Attempt to access to null thermo-table pointer");
            }
            return *thermoPtr_;
        }

        /*!\brief The transport-table.*/
        hur_nodiscard inline transportList &transTable() {
            if (!transportPtr_) {
                checkWarning("Attempt to access to null transport-table pointer");
            }
            return *transportPtr_;
        }

        /*!\brief The transport-table.*/
        hur_nodiscard inline const transportList &transTable() const {
            if (!transportPtr_) {
                checkWarning("Attempt to access to null transport-table pointer");
            }
            return *transportPtr_;
        }

        /*!\brief The reaction table.*/
        hur_nodiscard inline reactionList &reactions() {
            if (!reactionPtr_) {
                checkWarning("Attempt to access to null reaction pointer");
            }
            return *reactionPtr_;
        }

        /*!\brief The reaction table.*/
        hur_nodiscard inline const reactionList &reactions() const {
            if (!reactionPtr_) {
                checkWarning("Attempt to access to null reaction pointer");
            }
            return *reactionPtr_;
        }

#ifdef CUDA_PARALLEL

        /**
         * \brief Reaction table CUDA interface.
         */
        hur_nodiscard inline cuChem::cuChemInterface &reactionCUDA() {
            if (!cuChemInterfacePtr_) {
                checkWarning("Attempt to access to null reaction CUDA interface pointer");
            }
            return *cuChemInterfacePtr_;
        }

        /**
         * \brief Reaction table CUDA interface.
         */
        hur_nodiscard inline const cuChem::cuChemInterface &reactionCUDA() const {
            if (!cuChemInterfacePtr_) {
                checkWarning("Attempt to access to null reaction CUDA interface pointer");
            }
            return *cuChemInterfacePtr_;
        }

        /**
         * \brief The pointer of transport table in CUDA platform.
         */
        inline CUDATransport::transportListCUDA &transportCUDA() {
            if (!transportCUDAPtr_) {
                makeTransportTableCUDA();
            }
            return *transportCUDAPtr_;
        }

        /**
         * \brief The pointer of transport table in CUDA platform.
         */
        inline const CUDATransport::transportListCUDA &transportCUDA() const {
            if (!transportCUDAPtr_) {
                makeTransportTableCUDA();
            }
            return *transportCUDAPtr_;
        }

        /**
         * \brief The pointer of transport table in CUDA platform.
         */
        inline CUDATransport::transportListCUDA &transportCUDA(const cudaStreams &streams) {
            if (!transportCUDAPtr_) {
                makeTransportTableCUDA(streams);
            }
            return *transportCUDAPtr_;
        }

        /**
         * \brief The pointer of transport table in CUDA platform.
         */
        inline const CUDATransport::transportListCUDA &
        transportCUDA(const cudaStreams &streams) const {
            if (!transportCUDAPtr_) {
                makeTransportTableCUDA(streams);
            }
            return *transportCUDAPtr_;
        }

        inline void destroyTransportCUDA() const {
            if (transportCUDAPtr_) {
                transportCUDAPtr_->destroy();
            }
            transportCUDAPtr_.clear();
        }

        inline const cuChem::speciesTable &speciesCUDA() const {
            if (!speciesCUDAPtr_) {
                makeSpeciesCUDA();
            }
            return *speciesCUDAPtr_;
        }

        inline const cuChem::speciesTable &speciesCUDA(const cudaStreams &streams) const {
            if (!speciesCUDAPtr_) {
                makeSpeciesCUDA(streams);
            }
            return *speciesCUDAPtr_;
        }

        inline void destroySpeciesCUDA() const {
            if (speciesCUDAPtr_) {
                speciesCUDAPtr_->destroy();
            }
            speciesCUDAPtr_.clear();
        }

#endif // CUDA_PARALLEL

        /*!\brief Reutrn true if there is no reaction in the mixture.*/
        hur_nodiscard inline bool noReaction() const noexcept { return noReaction_; }

        /*!\brief Mixture average molecular weight field.*/
        hur_nodiscard cellRealArray W() const;

        /*!\brief Mixture average molecular weight.*/
        hur_nodiscard real W(const realArray &yi) const;

        /*!
         * \brief The density field for species i.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard cellRealArray rhoi(const cellRealArray &pi, const cellRealArray &Ti,
                                         const integer i, const bool onlyInternal = false) const;

        /*!
         * \brief The density field for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard cellRealArray rho(const cellRealArray &p, const cellRealArray &T,
                                        const bool onlyInternal = false) const;

        /**
         * \brief The density field for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard cellRealArray rho(const cellRealArray &p, const cellRealArray &T,
                                        const PtrList<cellRealArray> &_Yi,
                                        const bool onlyInternal = false) const;

        /**
         * \brief The density field for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        void rho(cellRealArray &rhom, const cellRealArray &p, const cellRealArray &T,
                 const PtrList<cellRealArray> &_Yi, const bool onlyInternal = false) const;

        /**
         * \brief The density field for species i.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard real rhoi(const real pi, const real Ti, const integer i) const;

        /**
         * \brief The density field for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard real rho(const real p, const real T, const realArray &yi) const;

        /**
         * \brief The density field for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard real rho(const real p, const real T, const integer cellI) const;

        /*!
         * \brief molecular viscosity for species i.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard cellRealArray mui(const cellRealArray &pi, const cellRealArray &Ti,
                                        const integer i, const bool onlyInternal = false) const;

        /*!
         * \brief molecular viscosity for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        cellRealArray mu(const cellRealArray &p, const cellRealArray &T,
                         const bool onlyInternal = false) const;

        /*
         * \brief Calculate the molecular viscosity for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[in] faceZoneId - Index of face zone in which the value should be calcualte.
         * \return Return the molecular viscosity for mixtures.
         */
        hur_nodiscard realArray mu(const cellRealArray &p, const cellRealArray &T,
                                   const integer faceZoneId) const;

        /*
         * \brief Calculate the molecular viscosity for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[out] mum - the molecular viscosity for mixtures
         * \param[in] onlyInternal - If true, only calculate the internal field.
         * \param[in] isExtrapolate - If true, the values in ghost cells are extrapolated from internal field and boundary field. Only effective when parameter onlyInternal is true.
         */
        void mu(const cellRealArray &p, const cellRealArray &T, cellRealArray &mum,
                const bool onlyInternal = false, const bool isExtrapolate = false) const;

        /*!
         * \brief Heat capacity at constant pressure for species i.
         * \param[in] rho - density [kg/m^3]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard cellRealArray cpi(const cellRealArray &rhoi, const cellRealArray &Ti,
                                        const integer i, const bool onlyInternal = false) const;

        /*!
         * \brief Heat capacity at constant pressure for mixtures.
         * \param[in] rho - density [kg/m^3]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard cellRealArray cp(const cellRealArray &rho, const cellRealArray &T,
                                       const bool onlyInternal = false) const;

        /*
         * \brief Heat capacity at constant pressure for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[in] faceZoneId - Index of face zone in which the value should be calcualte.
         * \return Return the heat capacity at constant pressure for mixtures.
         */
        hur_nodiscard realArray cp(const cellRealArray &p, const cellRealArray &T,
                                   const integer faceZoneId) const;

        /*
         * \brief Heat capacity at constant pressure for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[out] cpm -  The heat capacity at constant pressure for mixtures
         * \param[in] onlyInternal - If true, only calculate the internal field.
         * \param[in] isExtrapolate - If true, the values in ghost cells are extrapolated from internal field and boundary field. Only effective when parameter onlyInternal is true.
         */
        void cp(const cellRealArray &p, const cellRealArray &T, cellRealArray &cpm,
                const bool onlyInternal = false, const bool isExtrapolate = false) const;

        /*
         * \brief Enthalpy for species i.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[in] faceZoneId - Index of face zone in which the value should be calcualte.
         * \return Return the heat capacity at constant pressure for mixtures.
         */
        hur_nodiscard realArray hai(const cellRealArray &p, const cellRealArray &T,
                                    const integer faceZoneId, const integer isp) const;

        /*!
         * \brief Enthalpy for species i.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard cellRealArray hai(const cellRealArray &p, const cellRealArray &T,
                                        const integer i, const bool onlyInternal = false,
                                        const bool isExtrapolate = false) const;

        /*!
         * \brief Enthalpy for species i.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        void hai(const cellRealArray &p, const cellRealArray &T, cellRealArray &hi, const integer i,
                 const bool onlyInternal = false, const bool isExtrapolate = false) const;

        /*!
         * \brief Update absolute enthalpy for species i.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        void updateHai(const cellRealArray &p, const cellRealArray &T,
                       const bool onlyInternal = false, const bool isExtrapolate = false);

        /*!
         *\brief Ratio of specific heat capacities.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard real gamma(const real p, const real T, const realArray &yi) const;

        /*!
         * \brief Ratio of specific heat capacities.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard real gamma(const real p, const real T, const integer cellI) const;

        /*
         * \brief Ratio of specific heat capacities.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[in] faceZoneId - Index of face zone in which the value should be calcualte.
         * \return Return the ratio of specific heat capacities.
         */
        hur_nodiscard realArray gamma(const cellRealArray &p, const cellRealArray &T,
                                      const integer faceZoneId) const;

        /*
         * \brief Ratio of specific heat capacities.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[out] gm -  The ratio of specific heat capacities
         * \param[in] onlyInternal - If true, only calculate the internal field.
         * \param[in] isExtrapolate - If true, the values in ghost cells are extrapolated from internal field and boundary field. Only effective when parameter onlyInternal is true.
         */
        void gamma(const cellRealArray &p, const cellRealArray &T, cellRealArray &gm,
                   const bool onlyInternal = false, const bool isExtrapolate = false) const;

        //20210409 ��˼�� ����������gama��⺯��
        /*
         * \brief Ratio of specific heat capacities.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[out] gm -  The ratio of specific heat capacities
         * \param[in] isExtrapolate - If true, the values in ghost cells are extrapolated from internal field and boundary field. Only effective when parameter onlyInternal is true.
         */
        void gammaBoundary(const cellRealArray &p, const cellRealArray &T, cellRealArray &gm,
                           const bool isExtrapolate = false) const;

        /*!\brief Ratio of specific heat capacities.*/
        hur_nodiscard real gammai(const real p, const real T, const integer i) const;

        /*!\brief Ratio of specific heat capacities for ascoutic speed calculation.*/
        hur_nodiscard real gammaTh(const real p, const real rho, const real T,
                                   const realArray &yi) const;

        /*!\brief Ratio of specific heat capacities for ascoutic speed calculation.*/
        hur_nodiscard real gammaTh(const real p, const real rho, const real T,
                                   const integer cellI) const;

        /*!\brief Ratio of specific heat capacities for ascoutic speed calculation.*/
        hur_nodiscard real gammaThi(const real p, const real rho, const real T,
                                    const integer i) const;

        /*!
         * \brief Species diffusion for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         */
        hur_nodiscard PtrList<cellRealArray> Diff(const cellRealArray &p, const cellRealArray &T,
                                                  const bool onlyInternal = false) const;

        /*
         * \brief Calculate species diffusion for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[in] faceZoneId - Index of face zone in which the value should be calcualte.
         * \return Return species diffusion for mixtures.
         */
        hur_nodiscard realArrayArray Diff(const cellRealArray &p, const cellRealArray &T,
                                          const integer faceZoneId) const;

        /*
         * \brief Species diffusion for mixtures.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[out] Di -  Species diffusion for mixtures
         * \param[in] onlyInternal - If true, only calculate the internal field.
         * \param[in] isExtrapolate - If true, the values in ghost cells are extrapolated from internal field and boundary field. Only effective when parameter onlyInternal is true.
         */
        void Diff(const cellRealArray &p, const cellRealArray &T, PtrList<cellRealArray> &Di,
                  const bool onlyInternal = false, const bool isExtrapolate = false) const;

        /*
         * \brief Species diffusion for mixtures in ghost cells.
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[out] Di -  Species diffusion for mixtures
         * \param[in] isExtrapolate - If true, the values in ghost cells are extrapolated from internal field and boundary field.
         */
        void DiffBoundary(const cellRealArray &p, const cellRealArray &T,
                          PtrList<cellRealArray> &Di, const bool isExtrapolate = true) const;

    protected:
        void Diffi(realArray &Df, const cellRealArray &p, const cellRealArray &T,
                   const PtrList<cellRealArray> &yi, const integer cellI) const;

        hur_nodiscard inline realArray Diffi(const cellRealArray &p, const cellRealArray &T,
                                             const PtrList<cellRealArray> &yi,
                                             const integer cellI) const {
            realArray Df(species().size());
            Diffi(Df, p, T, yi, cellI);
            return Df;
        }

        void Diffi(realArray &Df, const real p, const real T, const realArray &yif) const;

    public:
        /*
         * \brief Calculate the molecular viscosity for mixtures and conefficients of thermal conductivity for mixtures.
         * \param[in] rho - density [kg/m^3]
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[out] mum - the molecular viscosity for mixtures
         * \param[out] kappam - the conefficients of thermal conductivity for mixtures
         * \param[in] onlyInternal - If true, only calculate the internal field.
         * \param[in] isExtrapolate - If true, the values in ghost cells are extrapolated from internal field and boundary field. Only effective when parameter onlyInternal is true.
         */
        void muKappa(const cellRealArray &rho, const cellRealArray &p, const cellRealArray &T,
                     cellRealArray &mum, cellRealArray &kappam, const bool onlyInternal = false,
                     const bool isExtrapolate = false) const;

        /*
         * \brief Calculate the molecular viscosity for mixtures and conefficients of thermal conductivity for mixtures in ghost cells.
         * \param[in] rho - density [kg/m^3]
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[out] mum - the molecular viscosity for mixtures
         * \param[out] kappam - the conefficients of thermal conductivity for mixtures
         * \param[in] isExtrapolate - If true, the values in ghost cells are extrapolated from internal field and boundary field.
         */
        void muKappaBoundary(const cellRealArray &rho, const cellRealArray &p,
                             const cellRealArray &T, cellRealArray &mum, cellRealArray &kappam,
                             const bool isExtrapolate = true) const;

        /*
         * \brief Calculate the molecular viscosity for mixtures and conefficients of thermal conductivity for mixtures.
         * \param[in] rho - density (dimensionless)
         * \param[in] p - static pressure (dimensionless)
         * \param[in] T - static temperature (dimensionless)
         * \param[out] mum - the molecular viscosity for mixtures (dimensionless)
         * \param[out] kappam - the conefficients of thermal conductivity for mixtures (dimensionless)
         * \param[in] faceZoneId - Index of face zone in which the value should be calcualte..
         */
        void muKappa(const cellRealArray &rho, const cellRealArray &p, const cellRealArray &T,
                     realArray &mum, realArray &kappam, const integer faceZoneId) const;

        /*
         * \brief Calculate the molecular viscosity for mixtures and conefficients of thermal conductivity for mixtures as well as the heat capacity at constant pressure for mixtures.
         * \param[in] rho - density [kg/m^3]
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[out] mum - the molecular viscosity for mixtures [kg/(m s)]
         * \param[out] kappam - the conefficients of thermal conductivity for mixtures [W / mK]
         * \param[out] cpm - the heat capacity at constant pressure for mixtures [J/(kg K)]
         * \param[in] onlyInternal - If true, only calculate the internal field.
         * \param[in] isExtrapolate - If true, the values in ghost cells are extrapolated from internal field and boundary field. Only effective when parameter onlyInternal is true.
         */
        void muKappaCp(const cellRealArray &rho, const cellRealArray &p, const cellRealArray &T,
                       cellRealArray &mum, cellRealArray &kappam, cellRealArray &cpm,
                       const bool onlyInternal = false, const bool isExtrapolate = false) const;

        /*
         * \brief Calculate the molecular viscosity for mixtures and conefficients of thermal conductivity for mixtures as well as the heat capacity at constant pressure for mixtures.
         * \param[in] rho - density [kg/m^3]
         * \param[in] p - static pressure [Pa]
         * \param[in] T - static temperature [K]
         * \param[out] mum - the molecular viscosity for mixtures [kg/(m s)]
         * \param[out] kappam - the conefficients of thermal conductivity for mixtures [W / mK]
         * \param[out] cpm - the heat capacity at constant pressure for mixtures [J/(kg K)]
         * \param[in] faceZoneId - Index of face zone in which the value should be calcualte..
         */
        void muKappaCp(const cellRealArray &rho, const cellRealArray &p, const cellRealArray &T,
                       realArray &mum, realArray &kappam, realArray &cpm,
                       const integer faceZoneId) const;

        void gasProperties(const cellRealArray &rho, const cellRealArray &p, const cellRealArray &T,
                           cellRealArray &mum, cellRealArray &kappam, cellRealArray &cpm,
                           PtrList<cellRealArray> &Di, PtrList<cellRealArray> &hai,
                           const bool onlyInternal = false, const bool isExtrapolate = false) const;

    protected:
        void gasProperties(const cellRealArray &rho, const cellRealArray &p, const cellRealArray &T,
                           realArray &mum, realArray &kappam, realArray &cpm, realArrayArray &Dif,
                           realArrayArray &haif, const integer faceZoneId) const;

    public:
        /*
         * \brief Calculate the molecular viscosity for mixtures.
         * \param[in] p - static pressure (dimensionless)
         * \param[in] T - static temperature (dimensionless)
         * \param[in] rho - density (dimensionless)
         * \param[in] v - velocity (dimensionless)
         * \param[out] E - total energy (dimensionless)
         * \param[in] internalOrBoundary - If true, only calculate the internal field, or only calculate the ghost field.
         * \param[in] isExtrapolate - If true, the values in ghost cells are extrapolated from internal field and boundary field. Only effective when parameter onlyInternal is true.
         */
        void E(const cellRealArray &p, const cellRealArray &T, const cellVectorArray &v,
               cellRealArray &E, const bool internalOrBoundary = false,
               const bool isExtrapolate = true) const;

        hur_nodiscard realArray E(const cellRealArray &p, const cellRealArray &T,
                                  const cellVectorArray &v, const integer faceZoneId) const;
    };
} // namespace OpenHurricane