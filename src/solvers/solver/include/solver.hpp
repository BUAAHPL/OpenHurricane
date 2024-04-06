/*!
 * \file solver.hpp
 * \brief Headers of base class of solver.
 *        The subroutines and functions are in the <i>solver.cpp</i> file.
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

#include "flowModel.hpp"
#include "monitors.hpp"
#include "patching.hpp"
#include "sourceTerms.hpp"
#include "spatialScheme.hpp"
#include "timeMarching.hpp"
#include "turbulenceModel.hpp"
#include "viscousFlux.hpp"
namespace OpenHurricane {
    // Forward declaration

    class writeFieldVars;
    class writeFaceZone;
    class calculateFieldVar;

    /*!\brief The base class of solver.*/
    class solver {
    public:
        /**\brief The bcType of time method.*/
        enum class timeTypes : short { UsualSolver, BDFUnsteadySolver };

    private:
        /*!\brief Require the auto patch?*/
        bool reqAutoPatch_;

        /*!\brief The auto patch?*/
        sharedPtrList<patching> autoPatch_;

        void setAutoPatch();

    protected:
        /**\brief The control of calculation.*/
        iteration &iter_;

        /**\brief The flow pointer.*/
        uniquePtr<flowModel> flowPtr_;

        /**\brief The spatial scheme pointer.*/
        uniquePtr<spatialScheme> invFluxPtr_;

        /**\brief Time method pointer.*/
        uniquePtr<timeMarching> timeMarcingPtr_;

        /**\brief Turbulence model pointer.*/
        uniquePtr<turbulenceModel> turbPtr_;

        /**\brief The timestep.*/
        //realArray dt_;
        cellRealArray dt_;

        /**\brief The number of cells in which the temperature beyond the lower limit.*/
        integer nLowTemp_;

        /**\brief The number of cells in which the temperature exceed the lower limit.*/
        integer nHighTemp_;

        /**\brief The number of cells in which the pressure beyond the lower limit.*/
        integer nLowPress_;

        /**\brief The number of cells in which the pressure exceed the lower limit.*/
        integer nHighPress_;

        bool hasTPLmtCells_;

    public:
        inline bool hasTPLmtCells() const noexcept;

        /**\brief The bcType of time method.*/
        timeTypes timeType_;

    public:
        inline void setUsualSolver() noexcept;
        inline void setBDFUnsteadySolver() noexcept;

    protected:
        /*!\brief Delete the pointer.*/
        template <class Type> inline void deletePointer(geometryArray<Type, cellMesh> *_gPtr) const;

        /*!\breif Linear exptrapolation variable to ghost field.*/
        template <class Type, class GeoMesh>
        void extrapolate(const Array<Array<Type>> &, geometryArray<Type, GeoMesh> &);

    public:
        // Calculate flux Jacobian

        /**
         * \brief The convective flux Jacobian.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[in] Vt - The contravariant velocity of the face of the control volume is set to zero for stationary grids
         * \param[out] AC - The convective flux Jacobian
         */
        virtual void Ac(const integer celli, const vector &normal, const real Vt,
                        realSquareMatrix &AC) const;

        /**
         * \brief The convective flux Jacobian.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[out] AC - The convective flux Jacobian
         */
        inline void Ac(const integer celli, const vector &normal, realSquareMatrix &AC) const;
        /**
         * \brief The derivative of pressure with respect to the conservative variables.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[out] pq - The derivative of pressure with respect to the conservative variables
         */
        virtual void dpdq(const integer celli, const vector &normal, realArray &pq) const;

        /**
         * \brief The convective flux Jacobian.
         * \param[in] celli - The index of cells
         * \param[in] normal - The unit normal vector of cell face
         * \param[in] dq - The increment of the conservative variables
         * \param[out] adq - The convective flux Jacobian multiplied by the increment of the conservative variables
         */
        virtual void Acdq(const integer celli, const vector &normal, const realArray &dq,
                          realArray &adq) const;

    protected:
        // Source terms

        uniquePtr<sourceTerms> sorcTermPtr_;

        void makeSorcTerm(const controller &sorcCont);

        inline sourceTerms &sorcTerm() noexcept;
        inline const sourceTerms &sorcTerm() const noexcept;

    public:
        declareClassName(solver);

        declareObjFty(solver, controller, (iteration & iter, const runtimeMesh &mesh),
                      (iter, mesh));

        /*
         * \brief Construct from components.
         * \param[in] iter - The iteration of computation.
         * \param[in] mesh - The computational mesh.
         */
        solver(iteration &iter, const runtimeMesh &mesh);

        hur_nodiscard static uniquePtr<solver> creator(iteration &iter, const runtimeMesh &mesh);

        /*!\brief Destructor.*/
        virtual ~solver() noexcept;

        /*!\brief Clear the solver.*/
        virtual void clear() noexcept;

        // Calculate

        /**\brief Solve the problem.*/
        virtual void solve();

        /**\brief Solve the problem using regular method.*/
        virtual void solving();

        /**\brief Solve the problem with BDF method.*/
        virtual void BDFSolve();

        /**\brief Update boundary conditions.*/
        virtual void bc();

        virtual void updateProperties() = 0;

        /**\brief Initialize the flow for the problem.*/
        inline virtual void initialize();

        /*!\brief To patch the given value for specific variable after initialization.*/
        bool initPatching();

        void autoPatching(const integer istep);

        /*!\brief Set RiemannValue::updated_ = false.*/
        inline virtual void iterRefresh();

        virtual void timeStep(realArray &dt);

        /**\brief Evaluate the convective flux.*/
        virtual void calculateFc();

        /**\brief Compute the viscous flux.*/
        virtual void calculateFv() = 0;

        /**\brief Caculate the source terms.*/
        virtual void calculateSource() = 0;

        /**\brief Other terms should be calculated.*/
        virtual void calculateOtherTerms() {}

        virtual void updatePrimitives(const bool shouldUpdateTemp = false);

        virtual void updateFlowOld();

        /**\brief Limit the temperature and pressure.*/
        virtual void limits();

        /**\brief Write the residuals and the results.*/
        virtual inline void write();

        /**\brief Conputational mesh.*/
        hur_nodiscard inline const runtimeMesh &mesh() const noexcept;

        /**\brief Iteration.*/
        hur_nodiscard inline const iteration &iter() const noexcept;

        /**\brief Time-step*/
        hur_nodiscard inline realArray &dt() noexcept;

        /**\brief Inviscous flux scheme.*/
        hur_nodiscard inline spatialScheme &invFlux() noexcept;

        /**\brief The time method.*/
        hur_nodiscard inline timeMarching &marching() noexcept;

        /**\brief The flag of shock.*/
        hur_nodiscard inline cellRealArray &shockFactor() noexcept;

        /**\brief The flag of shock.*/
        hur_nodiscard inline const cellRealArray &shockFactor() const noexcept;

        /**\brief Density*/
        hur_nodiscard inline cellRealArray &rho() noexcept;

        /**\brief Density*/
        hur_nodiscard inline const cellRealArray &rho() const noexcept;

        /**\brief Pressure*/
        hur_nodiscard inline cellRealArray &p() noexcept;

        /**\brief Pressure*/
        hur_nodiscard inline const cellRealArray &p() const noexcept;

        /**\brief Total energy.*/
        hur_nodiscard inline cellRealArray &E() noexcept;

        /**\brief Total energy.*/
        hur_nodiscard inline const cellRealArray &E() const noexcept;

        /**\brief Temperature.*/
        hur_nodiscard inline cellRealArray &T() noexcept;

        /**\brief Temperature.*/
        hur_nodiscard inline const cellRealArray &T() const noexcept;

        /**\brief Laminar dynamic viscosity.*/
        hur_nodiscard inline cellRealArray &mu() noexcept;

        /**\brief Laminar dynamic viscosity.*/
        hur_nodiscard inline const cellRealArray &mu() const noexcept;

        hur_nodiscard inline cellRealArray &mut() noexcept;

        hur_nodiscard inline const cellRealArray &mut() const noexcept;

        hur_nodiscard inline mixture &mixtures() noexcept;

        hur_nodiscard inline const mixture &mixtures() const noexcept;

        hur_nodiscard inline speciesList &specTable() noexcept;

        hur_nodiscard inline const speciesList &specTable() const noexcept;

        hur_nodiscard inline cellVectorArray &v() noexcept;

        hur_nodiscard inline const cellVectorArray &v() const noexcept;

        /**\brief Laminar Prandtl number.*/
        hur_nodiscard inline real prl() const noexcept;

        /**\brief Turbulent Prandtl number*/
        hur_nodiscard inline real prt() const noexcept;

        /*!\brief Schmidt number.*/
        hur_nodiscard inline real sct() const noexcept;

        hur_nodiscard inline cellRealArray &kappal() noexcept;

        hur_nodiscard inline const cellRealArray &kappal() const noexcept;

        hur_nodiscard inline cellRealArray &kappat() noexcept;

        hur_nodiscard inline const cellRealArray &kappat() const noexcept;

        hur_nodiscard inline cellRealArray kappaEff() const;

        hur_nodiscard inline cellRealArray &cp() noexcept;

        hur_nodiscard inline const cellRealArray &cp() const noexcept;

        hur_nodiscard inline cellRealArray &gama() noexcept;

        hur_nodiscard inline const cellRealArray &gama() const noexcept;

        /*!\brief The mass fractions of species.*/
        hur_nodiscard inline PtrList<cellRealArray> &yi() noexcept;

        /*!\brief The mass fractions of species.*/
        hur_nodiscard inline const PtrList<cellRealArray> &yi() const noexcept;

        /*!\brief The absolute enthalpy of species.*/
        hur_nodiscard inline PtrList<cellRealArray> &hi() noexcept;

        /*!\brief The absolute enthalpy of species.*/
        hur_nodiscard inline const PtrList<cellRealArray> &hi() const noexcept;

        /*!\brief The mass diffusion coefficient of species.*/
        hur_nodiscard inline PtrList<cellRealArray> &Diff() noexcept;

        /*!\brief The mass diffusion coefficient of species.*/
        hur_nodiscard inline const PtrList<cellRealArray> &Diff() const noexcept;

        /**\brief The flag of state to evaluate temperature.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        hur_nodiscard inline integerArray &temperatureFlag() noexcept;

        /**\brief The flag of state to evaluate temperature.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        hur_nodiscard inline const integerArray &temperatureFlag() const noexcept;

        /**\brief The flag of state to evaluate pressure.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        hur_nodiscard inline const integerArray &pressureFlag() const noexcept;

        /**\brief The flag of state to evaluate pressure.
         * -# 1 -  Regular state.
         * -# 0 -  Beyond the lower limit.
         * -# 2 -  Exceed the higher limit.
         */
        hur_nodiscard inline integerArray &pressureFlag() noexcept;

        /*!\brief The flag to adjust CFL number.
         * -# 1 - Do not adjust CFL number.
         * -# 0 - Reduce the CFL number.
         */
        hur_nodiscard inline const cellIntegerArray &CFLFlag() const noexcept;

        /*!\brief The flag to adjust CFL number.
         * -# 1 - Do not adjust CFL number.
         * -# 0 - Reduce the CFL number.
         */
        hur_nodiscard inline cellIntegerArray &CFLFlag() noexcept;

        inline virtual void firstSolveSplittingSource(const real phyDt);
        inline virtual void secondSolveSplittingSource(const real phyDt);

        inline virtual void previousSource();
        inline virtual void postSource();

    public:
        template <class Type>
        hur_nodiscard Array<Type> gatherCellDataInMaster(const Array<Type> &cellQ) const;

    protected:
        bool useLowMachPrecon_;

    public:
        hur_nodiscard inline bool useLowMachPrecon() const noexcept;

    protected:
        // Low mach preconditioning parameters
        mutable cellRealArray *a4Ptr_;

    public:
        inline cellRealArray &a4();
        inline const cellRealArray &a4() const;

    protected:
        mutable cellRealArray *a5Ptr_;

    public:
        inline cellRealArray &a5();
        inline const cellRealArray &a5() const;

    protected:
        /**
         * \brief The transformation matrix from the primitive into the conservative variables.
         * \param[in] cellI - The index of cell.
         * \param[out] PM - The transformation matrix
         */
        virtual void PP(const integer cellI, realSquareMatrix &PM1, realSquareMatrix &PP,
                        realSquareMatrix &FM1, realSquareMatrix &FP, real &aa4, real &aa5) const;
    };
} // namespace OpenHurricane

#include "solver.inl"