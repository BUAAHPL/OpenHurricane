/*!
 * \file chemistrySource.hpp
 * \brief Header of chemistry source.
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

#include "ODEsSolver.hpp"
#include "flowModel.hpp"

#ifdef CUDA_PARALLEL
#include "cudaStreams.hpp"
#include "nThreadsAndBlocks.hpp"
#endif // CUDA_PARALLEL

namespace OpenHurricane {
    /**\brief The class of chemistry source.*/
    class chemistrySource {
    public:
        enum class ReactionFlowType : short { ConstantPressure, ConstantVolume, IsolatedSystem };

    protected:
        /*!\brief Flows.*/
        flowModel &flows_;

        /*!\brief Reference to the runtime mesh.*/
        const runtimeMesh &mesh_;

        /*!\brief Reference to the field of species mass fractions.*/
        PtrList<cellRealArray> &yi_;

        /*!\brief Reaction table.*/
        reactionList &reactions_;

        /*!\brief Species table.*/
        const speciesList &species_;

        /*!\brief the number of species.*/
        const integer nsp_;

        /*!\brief The number of activate species of current mechanism.*/
        integer nSpc_;

        /*!\brief The number of reactions.*/
        const integer nrc_;

        /*!\brief The number of activate reactions of current mechanism.*/
        integer nRact_;

        /*!\brief Is calculating chemistry kinetics?*/
        bool isReacting_;

        /** \brief Is using simplified Jacobians? */
        bool simplifiedJacobians_;

        /*!\brief Temporary concentration array.*/
        mutable realArray ci_;
        mutable realArray yiTmp_;

        /*!\brief Temporary rate-of-change of concentration field.*/
        mutable realArray dcdt_;

        /*!\brief Chemical source terms per specie [kg/m^3/s].*/
        mutable realArrayArray *RiPtr_;

    public:
        /*!\brief Chemical source terms [kg/m^3/s].*/
        inline realArrayArray &Ri();

        /*!\brief Chemical source terms per specie [kg/m^3/s].*/
        inline const realArrayArray &Ri() const;
#ifdef TEST_PROCESS_TIME
    protected:
        integer odeCountIter_;

        virtual inline void setODECountIter() noexcept;

    public:
        hur_nodiscard inline integer odeCountIter() const noexcept;

#endif
    protected:
        /*!\brief Temporary rate-of-change of concentration field.*/
        mutable realArray *dtInitPtr_;

        realArray &dtInit(const realArray &dt, const real dtFactor = real(1));

        real dtMax_;

        /**\brief ha0 or e0.*/
        real hea0_;

        real T_;

        real p0_;
        real p_;

        real rho0_;
        real rho_;

        ReactionFlowType reactionFlowType_;

#ifdef TEST_PROCESS_TIME

        /**
         * \brief The time consumed by solving chemical ODEs.
         */
        real calcTime_;

        /**
         * \brief The time cost by reducing chemical kinetics.
         */
        real reducionTime_;

    public:
        /**
         * \brief The time consumed by solving chemical ODEs.
         */
        hur_nodiscard inline real totalSolveIngTime() const noexcept;

        /**
         * \brief The time cost by reducing chemical kinetics.
         */
        hur_nodiscard inline real reducionTime() const noexcept;

        inline virtual void writeTimeTitle(fileOsstream &fos) const;
        inline virtual void writeSolveODETime(fileOsstream &fos) const;

#endif // TEST_PROCESS_TIME

    protected:
        /**
         * !\brief The number of spcies in every cell.
         */
        mutable cellRealArray *nSPInEveryCellPtr_;

        /**
         * !\brief The number of reactions in every cell.
         */
        mutable cellRealArray *nRCInEveryCellPtr_;

    public:
declareClassNames;
declareObjFty(chemistrySource,controller,
                            (flowModel & flows, const controller &cont), (flows, cont));

        // Constructors

        /*!\brief Disallow null constructor.*/
        chemistrySource() = delete;

        /**
         * \brief Disallow bitwise copy constructor.
         */
        chemistrySource(const chemistrySource &) = delete;

        /*!\brief Construct from flow and controller.*/
        chemistrySource(flowModel &flows, const controller &cont);

        static uniquePtr<chemistrySource> creator(flowModel &flows, const controller &cont);

        /*!\brief Destructor.*/
        virtual ~chemistrySource() noexcept;

        /*!\brief Reference to the runtime mesh.*/
        inline const runtimeMesh &mesh() const noexcept;

        /*!\brief Reference to the field of species mass fractions.*/
        inline PtrList<cellRealArray> &yi() noexcept;

        /*!\brief Reference to the field of species mass fractions.*/
        inline const PtrList<cellRealArray> &yi() const noexcept;

        /*!\brief Reaction table.*/
        inline reactionList &reactions() noexcept;

        /*!\brief Reaction table.*/
        inline const reactionList &reactions() const noexcept;

        /*!\brief Species table.*/
        inline const speciesList &species() const noexcept;

        /*!\brief The number of species.*/
        virtual inline integer nsp() const noexcept;

        /*!\brief The number of activate species of current mechanism.*/
        hur_nodiscard inline integer nSpc() const noexcept;

        /*!\brief The number of reactions.*/
        virtual inline integer nrc() const noexcept;

        /*!\brief The number of activate reactions of current mechanism.*/
        hur_nodiscard inline integer nRact() const noexcept;

        /*!\brief Is calculating chemistry kinetics?*/
        inline bool isReacting() const noexcept;

        /*!\brief The number of equations.*/
        inline virtual integer nEqns() const;
        inline virtual integer kjkj(integer j) { return j; }

        inline void setReactionFlowType(const ReactionFlowType t) noexcept;

        inline void setConstantPressure() noexcept;

        inline void setConstantVolume() noexcept;

        inline bool isConstPressure() const noexcept;
        inline bool isConstVolume() const noexcept;

        inline void init(const real hOre, const real T0, const real p0, const real rho0) noexcept;

        /*!\brief Return the absolute enthalpy before reaction.*/
        inline real hea0() const noexcept;

        /*!\brief Return the absolute enthalpy before reaction.*/
        inline real &hea0() noexcept;

        /*!\brief Return the pressure before reaction.*/
        inline real p0() const noexcept;

        /*!\brief Return the pressure before reaction.*/
        inline real &p0() noexcept;

        /*!\brief Return the temperature after reaction.*/
        inline real T() const noexcept;

        /*!\brief Return the temperature after reaction.*/
        inline real &T() noexcept;

        /*!\brief Return the pressure after reaction.*/
        inline real p() const noexcept;

        /*!\brief Return the pressure after reaction.*/
        inline real &p() noexcept;

        /*!\brief Return the density before reaction.*/
        inline real rho0() const noexcept;

        /*!\brief Return the density before reaction.*/
        inline real &rho0() noexcept;

        /*!\brief Return the density after reaction.*/
        inline real rho() const noexcept;

        /*!\brief Return the density after reaction.*/
        inline real &rho() noexcept;

        /*!\brief Return the size of the species list.*/
        /*inline integer speciesSize() const noexcept;*/

        /*!\brief The mixture.*/
        inline const mixture &mixtures() const;

        /*!\brief The mixture.*/
        inline mixture &mixtures();

        /*!\brief The thermo-table.*/
        inline const thermoList &therm() const;

        /*!\brief The thermo-table.*/
        inline thermoList &therm();

        /*!\brief Return the net rate of j_th reaction.
         * \param[in] p - static pressure [Pa].
         * \param[in] T - static temperature [K].
         * \param[in] c - molar concentrations of species [kmol/m^3].
         * \param[in] j - the index of reactions.
         * \return qj - the net rate of j_th reaction.
         */
        //real qj(const real p, const real T, const realArray& c, const integer j);

        /*!\brief Return the net rate of j_th reaction with the diagonal Jacobian matrix.
         * \param[in] p - static pressure [Pa].
         * \param[in] T - static temperature [K].
         * \param[in] c - molar concentrations of species [kmol/m^3].
         * \param[in] j - the index of reactions.
         * \param[out] diagJac - the diagonal Jacobian matrix for the first Ns - 1 species.
         * \return qj - the net rate of j_th reaction.
         */
        real qj(const real p, const real T, const realArray &c, const integer j, realArray &diagJac,
                const bool calcLastSpeciesJac = false);

        /*!\brief Return the net rate of j_th reaction with the forward and reverse rate, respectively.
         * \param[in] p - static pressure [Pa].
         * \param[in] T - static temperature [K].
         * \param[in] c - molar concentrations of species [kmol/m^3].
         * \param[in] j - the index of reactions.
         * \param[out] qfj - forward rate of j_th reaction.
         * \param[out] qrj - reverse rate of j_th reaction.
         * \return qj - the net rate of j_th reaction.
         */
        /*real qfrj
        (
                const real p,
                const real T,
                const realArray& c,
                const integer j,
                real& qfj,
                real& qrj
        );*/

        /*!\brief Return the net rate of j_th reaction with the diagonal Jacobian matrix
         *        and the forward and reverse rate, respectively.
         * \param[in] p - static pressure [Pa].
         * \param[in] T - static temperature [K].
         * \param[in] c - molar concentrations of species [kmol/m^3].
         * \param[in] j - the index of reactions.
         * \param[out] qfj - forward rate of j_th reaction.
         * \param[out] qrj - reverse rate of j_th reaction.
         * \param[out] diagJac - the diagonal Jacobian matrix.
         * \param[in] calcLastSpeciesJac - If true, the size of diagJac is Ns, otherwise Ns-1.
         * \return qj - the net rate of j_th reaction.
         */
        real qfrj(const real p, const real T, const realArray &c, const integer j, real &qfj,
                  real &qrj, realArray &diagJac, const bool calcLastSpeciesJac = true);

        /*!\brief Chemical source term for chemical ODEs.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \return The chemical source term with an array size of Ns [in molar unit]
         */
        virtual realArray omega(const real p, const real T, const realArray &c);

        /*!\brief Chemical source term for chemical ODEs.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \param[out] w - the chemical source term with an array size of Ns [in molar unit]
         */
        virtual void omega(const real p, const real T, const realArray &c, realArray &w);

        /*!\brief Chemical source term for all species and Jacobian of source terms.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \param[out] w - the chemical source term with an array size of Ns [in molar unit]
         * \param[out] dwdy - {d_w/d_y: y = (c_1,...,c_Ns,T)} the chemical source term Jacobian with an array size of Ns+1 x Ns+1 [in time unite]
         */
        virtual void omegaFullJacobian(const real p, const real T, const realArray &c, realArray &w,
                                       realSquareMatrix &dwdy);

        /*!\brief Chemical source term for chemical ODEs.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \param[out] w - the chemical source term with an array size of Ns [in molar unit]
         */
        virtual void omega(const real rho, const real p, const real T, const realArray &c,
                           const realArray &yi, realArray &w, realArray &diagJac);

        /*!\brief Chemical source terms for coupled NS equations.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \param[out] w - the chemical source term with an array size of Ns - 1 [in molar unit]
         */
        virtual void omegaCoupled(const real p, const real T, const realArray &c, realArray &w);

        /*!\brief Chemical source terms for coupled NS equations.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \param[out] w - the chemical source term with an array size of Ns - 1 [in molar unit]
         * \param[out] diagJac - The Jacobian diagonal terms of an array size of Ns - 1
         */
        virtual void omegaCoupled(const real rho, const real p, const real T, const realArray &c,
                                  const realArray &yi, realArray &w, realArray &diagJac);

    protected:
        realArray *diagJacPerReacPtr_;

        inline realArray &diagJacPerReac();

    public:
        /*!\brief Use the approach based on the forward reaction rates.
         *		  Chemical source terms for coupled NS equations.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \param[out] omega - the chemical source term with an array size of Ns - 1 [in molar unit]
         * \return The characteristic chemical time scale [s]
         */
        virtual real omegaCoupledAndtc(const real p, const real T, const realArray &c,
                                       realArray &omega);

        /*!\brief Use the approach based on the forward reaction rates.
         *		  Chemical source terms for coupled NS equations.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \param[out] omega - the chemical source term with an array size of Ns - 1 [in molar unit]
         * \param[out] diagJac - The Jacobian diagonal terms of an array size of Ns - 1
         * \return The characteristic chemical time scale [s]
         */
        virtual real omegaCoupledAndtc(const real rho, const real p, const real T,
                                       const realArray &c, const realArray &yi, realArray &omega,
                                       realArray &diagJac);

        /*!\brief The time scale for all species using the approach based on the species formation rates.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \param[out] ts - the time scale for all speccies [s]
         */
        virtual void timescale(const real p, const real T, const realArray &c, realArray &ts);

        /*!\brief The characteristic chemical time scale for all species using the approach based on the species formation rates.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         */
        virtual hur_nodiscard real tcSFR(const real p, const real T, const realArray &c);

        /*!\brief The time scale for all species using the approach based on the Jacobian diagonal terms without the last species.
         * \param[in] rho - density [kg/m^3]
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \param[in] yyi - species mass fraction
         * \param[out] ts - the time scale for all speccies except the last species [s]
         */
        virtual void timescale2(const real rho, const real p, const real T, const realArray &c,
                                const realArray &yyi, realArray &ts);

        /*!\brief The time scale for all species using the approach based on the Jacobian diagonal terms.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3] (Should be full size of all species)
         * \param[out] ts - the time scale for all speccies [s]
         */
        virtual void timescale2(const real p, const real T, const realArray &c, realArray &ts);

        /*!\brief The time scale for all reaction modes using the approach based on the Jacobian eigen decomposition.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3] (Should be full size of all species)
         * \param[out] ts - the time scale for all reaction modes [s]
         */
        virtual void timescaleReacMode(const real p, const real T, const realArray &c,
                                       realArray &ts);

        /*!\brief The characteristic chemical time scale using the approach based on the forward reaction rates.
         * \param[in] p - pressure [Pa]
         * \param[in] T - temperature [K]
         * \param[in] c - molar concentration [kmol/m^3]
         * \param[out] ts - the time scale for all speccies [s]
         */
        virtual real tc(const real p, const real T, const realArray &c);

        /**\brief Calculate the derivatives in dydt.*/
        virtual void DyDt(const real t, const realArray &y, realArray &dydt);

    protected:
        /**\brief Calculate the derivatives in dydt.*/
        virtual void DyDtSlected(const real t, const realArray &y, realArray &dydt,
                                 const bool withTemp);

    public:
        virtual void updateOtherValue(const realArray &c);

        /*!\brief Make sure the summation of mass fractions of all species equal to 1.
         * \param[in,out] yi - mass fraction
         */
        virtual void correctY(realArray &y);

        /**
         *\brief Calculate the Jacobian of the system.
         *  Need by the stiff-system solvers.
         */
        virtual void jacobian(const real t, const realArray &ci, realArray &dcidt,
                              realSquareMatrix &dfdci);

    protected:
        /**
         *\brief Calculate the Jacobian of the system.
         *  Need by the stiff-system solvers.
         */
        virtual void jacobianSelected(const real t, const realArray &ci, realArray &dcidt,
                                      realSquareMatrix &dfdci, const bool withTemp);

    protected:
        /**
         * \brief Temperature positivity limit rate. Default is 0.2.
         */
        real TemperatureAlpha_;

        /**
         * \brief Temperature time step reduction. Default is 0.25.
         */
        real TemperatureBeta_;

        real TemperatureThreshold_;

    public:
        /**
         * \brief Check the new y is acceptable.
         * \param[in,out] dt - The timestep. If unacceptable, then reduce dt.
         * \param[in] yOld - The old solution array
         * \param[in] yNew - The new solution array
         * \return True if acceptable.
         */
        virtual bool checkNewY(real &dt, const realArray &yOld, const realArray &yNew);

        virtual void DomegaDciT(const real rho, const real T, const real p, const realArray &yi,
                                realArray &omegaj, realSquareMatrix &dwjdciT) const;

        /*!\brief Convert molar concentration to mass fraction
         * \param[in] ci - the molar concentration
         * \param[out] rho - the density
         * \return The array of mass fraction
         */
        realArray ciToYi(const realArray &ci, real &rho) const;

        void ciToYi(const realArray &ci, const real rho, realArray &yi) const;
        real ciToYi(const realArray &ci, realArray &yi) const;

        /*!\brief Convert molar concentration to mass fraction
         * \param[in] ci - the molar concentration
         * \return The array of mass fraction
         */
        realArray ciToYi(const realArray &ci) const;

        /*!\brief Convert mass fraction to molar concentration.
         * \param[in] yi - the mass fraction of species
         * \param[in] rho - the density
         * \return The array of molar concentration
         */
        realArray yiToCi(const realArray &yi, const real rho) const;

        /*!\brief Convert mass fraction to molar concentration.
         * \param[in] yi - the mass fraction of species
         * \param[in] rho - the density
         * \param[out] ci - the array of molar concentration
         */
        void yiToCi(const realArray &yi, const real rho, realArray &ci) const;

    protected:
        hur_nodiscard virtual real getTemperature(const realArray &ci,
                                                  const bool check = false) const;

    public:
        /*!\brief Solve chemical ODEs.
         * \param[in] t - the time for reacting
         * \param[in] dT - the sub-timestep for reaction
         * \param[in] p - the pressure
         * \param[in] T - the temperature
         * \param[in] rho - the density
         * \param[in,out] yi - the mass fraction of species
         * \return The last sub-timestep used.
         */
        real solve(const real t, const real dT, const real p, const real T, const real rho,
                   realArray &yi);

        virtual real solveTest(const real t, const real dT, real &_p, real &_T, real &_rho,
                               realArray &yi, fileOsstream &fos);

#ifdef TEST_PROCESS_TIME

        /**
         * \brief Only used for test solveTest().
         * \return 0 - for not found in table; 1 - for found in table
         */
        inline virtual hur_nodiscard integer findInTable() const noexcept { return 0; }
#endif
        /*!\brief Solve chemical ODEs.
         * \param[in,out] yi - the mass fraction of species
         * \param[in,out] dT - the time for reaction
         * \param[in,out] subDT - the sub-timestep
         */
        virtual void solve(realArray &yi, real &dT, real &subDT);

        virtual void calculateSourceTerms(const bool withoutLastSpec = true);
        virtual void calculateSourceTermsImp(const bool withoutLastSpec = true);

        /**
         * \brief Chemical source term for all species and Jacobian of source terms at cellI.
         * \param[in] cellI - Cell index
         * \param[out] dwdy - {d_w/d_y: y = (c_1,...,c_Ns,T)} the chemical source term Jacobian with an array size of Ns+1 x Ns+1 [in time unite]
         */
        virtual void sourceTermsForCellI(const integer cellI, realSquareMatrix &dwdy);

        /**
         * \brief Return heat release rate [J/(m^3 s)].
         */
        virtual realArray heatReleaseRate();

        /**
         * \brief Return heat release rate [J/(m^3 s)].
         */
        virtual real heatReleaseRate(const real p, const real T, const realArray &c);

        /*!\brief Solve chemical source terms.
         * \param[in] dt - timestep of every cell for reacting.
         * \param[in] dtFactor - the facetor to scale dt.
         * \return Return the minimum timestep.
         */
        virtual real solve(const realArray &dt, const real dtFactor);

        /*!\brief Solve chemical source terms and update species mass fraction and temperature.
         * \param[in] dt - timestep of every cell for reacting.
         * \param[in] dtFactor - the facetor to scale dt.
         */
        virtual void solveUpdate(const realArray &dt, const real dtFactor);

        /*virtual void solve(
                cellRealArray& rhoi,
                cellRealArray& pi,
                cellRealArray& Ti,
                PtrList<cellRealArray>& yii,
                const realArray& dt,
                realArray& subDt,
                const real dtFactor
        );*/

        /*!\brief Return access to the chemical source terms for species i. [kg/(m^3 s)]*/
        inline realArray &Ri(const integer i);

        /*!\brief Return const access to the chemical source terms for species i. [kg/(m^3 s)]*/
        inline const realArray &Ri(const integer i) const;

        inline virtual bool isReduced() const noexcept;
#ifdef CUDA_PARALLEL

        virtual void calculateSourceTermsAsync(real *hur_restrict hostYiRhoTPtr_,
                                               cu2DArray<cu_real> &dYiRhoT,
                                               real *hur_restrict RRi, cu2DArray<cu_real> &dRi,
                                               const cudaStreams &streams);

        virtual void calculateSourceTermsImpAsync(
            real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT,
            real *hur_restrict RRi, cu2DArray<cu_real> &dRi, real *hur_restrict dRidrhoyi,
            cu2DArray<cu_real> &ddRidrhoyi, const cudaStreams &streams);

        virtual void calculateSourceTermsImpAsyncHybrid(
            real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT,
            real *hur_restrict RRi, cu2DArray<cu_real> &dRi, cu_float *hur_restrict dRidrhoyi,
            cu2DArray<cu_float> &ddRidrhoyi, const cudaStreams &streams);

        
        virtual void calculateSourceTermsAsyncSlice(real *hur_restrict hostYiRhoTPtr_,
                                                    cu2DArray<cu_real> &dYiRhoT,
                                                    real *hur_restrict RRi,
                                                    cu2DArray<cu_real> &dRi,
                                                    const cudaStreams &streams, const integer pos,
                                                    const integer offset) {
            LFatal("This function is not implemented");
        }

        virtual void calculateSourceTermsImpAsyncSlice(
            real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT,
            real *hur_restrict RRi, cu2DArray<cu_real> &dRi, real *hur_restrict dRidrhoyi,
            cu2DArray<cu_real> &ddRidrhoyi, const cudaStreams &streams, const integer pos,
            const integer offset) {
            LFatal("This function is not implemented");
        }

        virtual void calculateSourceTermsImpAsyncHybridSlice(
            real *hur_restrict hostYiRhoTPtr_, cu2DArray<cu_real> &dYiRhoT,
            real *hur_restrict RRi, cu2DArray<cu_real> &dRi, cu_float *hur_restrict dRidrhoyi,
            cu2DArray<cu_float> &ddRidrhoyi, const cudaStreams &streams, const integer pos,
            const integer offset) {
            LFatal("This function is not implemented");
        }

        virtual void createReactionTable() const;
        virtual void createReactionTableAsync(const cudaStreams &streams) const;

        virtual void destroyReactionTable() const;

#endif // CUDA_PARALLEL

        /**
         * !\brief The number of spcies in every cell.
         */
        inline virtual const cellRealArray &nSPInEveryCell() const;

        /**
         * !\brief The number of spcies in every cell.
         */
        inline virtual cellRealArray &nSPInEveryCell();

        /**
         * !\brief The number of reactions in every cell.
         */
        inline virtual const cellRealArray &nRCInEveryCell() const;

        /**
         * !\brief The number of reactions in every cell.
         */
        inline virtual cellRealArray &nRCInEveryCell();

    public:
        virtual void setCellLoadWeights();

    protected:
        mutable cellRealArray *cellSourceCalTimePtr_;

    public:
        hur_nodiscard inline const cellRealArray &cellSourceCalTime() const;
        hur_nodiscard inline cellRealArray &cellSourceCalTime();

#ifdef TEST_PROCESS_TIME
        inline virtual hur_nodiscard integer nSpcGroup() const noexcept;
        inline virtual hur_nodiscard integerList nEachSpcGroup() const noexcept;

#endif // TEST_PROCESS_TIME
    };
} // namespace OpenHurricane

#include "chemistrySource.inl"