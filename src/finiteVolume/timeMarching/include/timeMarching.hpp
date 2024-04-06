/*!
 * \file timeMarching.hpp
 * \brief Headers of class of the time marching.
 *        The subroutines and functions are in the <i>timeMarching.cpp</i> file.
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
#include "fVArraysInclude.hpp"
#include "flowModel.hpp"
#include "matrices.hpp"
#include "pseudoTime.hpp"
#include "smartPointerList.hpp"
#include "spatialScheme.hpp"

namespace OpenHurricane {
    /* Forward declaration */

    class solver;

    /*!\brief The base class of time marching.*/
    class timeMarching {

    public:
        /**
         * \brief The method type for source terms treating.
         */
        enum class sourceTermTreating : short {
            EXPLICITSOURCE,
            KIMDIAGIMPLICITSOURCE,
            FULLJACOBIAN,
            FULLJACOBIANTABLE
        };

    private:
        /*!\brief Hold reference to the runtime mesh.*/
        const runtimeMesh &mesh_;

        /*!\brief The total number of the primitive parameters*/
        integer nPrimitives_;

    protected:
        /** \brief Reference to solver */
        solver &solver_;

        /*!\brief Velocity field.*/
        cellVectorArray &v_;

        /*!\brief Thermo.*/
        const flowModel &flow_;

        /*!\brief The list of primitive objects.*/
        List<object *> objectList_;

        /** \brief The number of scalar primitive parameters. */
        integer countParam_;

        integerVector2DList paramMap_;

        /*!\brief The increment.*/
        realArrayArray dq_;

        /**
         * \brief Transfer dq in cut face and periodic face boundaries.
         */
        inline void transferDq();

        /*!\brief Time step.*/
        realArray &dt_;

        /**\brief Timestep reduced factor.*/
        real dtReduceFactor_;

        cellRealArray &shockFactor_;

        bool isEuler_;

        /*!\brief The spectral radius of the convective fluxes*/
        faceVector2DArray rai_;

        /*!\brief The spectral radius of the viscous fluxes.*/
        faceVector2DArray rav_;

        /** \brief Pseudo time step */
        uniquePtr<pseudoTime> pseudoTimePtr_;

        void calcFluxSpectRadius();

        /** \brief Full point-implicit Jacobian matrix array pointer. */
        uniquePtr<cellRealSquareMatrixArray> JacPtr_;

    protected:
        /** \brief The method of source terms treating. */
        sourceTermTreating sourceTermTreat_;

    public:
        declareClassName(timeMarching);
        declareObjFty(timeMarching, controller,
                      (const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
                       solver &_solver, cellVectorArray &v),
                      (cont, mesh, flowM, _solver, v));

        /*!\brief Disallow null constructor.*/
        timeMarching() = delete;

        /*!\brief Disallow copy constructor.*/
        timeMarching(const timeMarching &) = delete;

        timeMarching(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
                     solver &_solver, cellVectorArray &v);

        static uniquePtr<timeMarching> creator(const controller &cont, const runtimeMesh &mesh,
                                               const flowModel &flowM, solver &_solver,
                                               cellVectorArray &v);

        /*!\brief Destructor.*/
        virtual ~timeMarching() noexcept;

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept;

        hur_nodiscard const iteration &iter() const noexcept;

        /*!\brief Time step.*/
        hur_nodiscard inline realArray &dt() noexcept;

        /*!\brief The spectral radius of the convective fluxes*/
        hur_nodiscard inline const faceVector2DArray &rai() const noexcept;

        /*!\brief The spectral radius of the viscous fluxes.*/
        hur_nodiscard inline const faceVector2DArray &rav() const noexcept;

        hur_nodiscard inline bool isEuler() const noexcept;

        /**
         * \brief If use explicit source?.
         */
        hur_nodiscard inline bool explicitSource() const noexcept;

        hur_nodiscard inline bool diagonalImpSource() const noexcept;

        hur_nodiscard inline bool fullJacobianSource() const noexcept;

        hur_nodiscard inline bool fullJacobianSourceTable() const noexcept;

        /** \brief Full point-implicit Jacobian matrix array. */
        hur_nodiscard inline cellRealSquareMatrixArray &Jacobian();

        /** \brief Full point-implicit Jacobian matrix array. */
        hur_nodiscard inline const cellRealSquareMatrixArray &Jacobian() const;

        hur_nodiscard inline const cellRealArray &shockFactor() const noexcept;
        hur_nodiscard inline cellRealArray &shockFactor() noexcept;

        /*!\brief The total number of the primitive parameters*/
        hur_nodiscard inline integer nPrimitives() const noexcept;

        /** \brief The number of scalar primitive parameters. */
        hur_nodiscard inline integer countParam() const noexcept;

        hur_nodiscard inline realArrayArray &dq() noexcept;

        inline const hur_nodiscard realArrayArray &dq() const noexcept;

        template <class Type> void addObject(geometryArray<Type, cellMesh> &);

        hur_nodiscard real vistParam(const integer paramI, const integer cellI) const;

        hur_nodiscard realArray vistParam(const integer cellI) const;

        void pushParam(const integer paramI, const integer cellI, const real value) const;

        /**\brief Get access to the right hand side of equation paramI by giving the index of the cell.*/
        hur_nodiscard real vistRhs(const integer paramI, const integer cellI) const;

        /**\brief Get access the right hand side of all equations by giving the cell index.*/
        hur_nodiscard realArray vistRhs(const integer cellI) const;

        /**\brief Get access to the diagonal of Jacobian matrices by giving the index of the parameter and the cell.*/
        hur_nodiscard real vistDiagSource(const integer paramI, const integer cellI) const;

        /**\brief Get access to the diagonal of Jacobian matrices by giving the cell index.*/
        hur_nodiscard realArray vistDiagSource(const integer cellI) const;

        virtual void initializing();

        /**\brief Store the last tiem step primitive flow variables.*/
        void storeOldPrimitiveValue();

        /**\brief To compute the time-step for time marching.*/
        virtual void timeStep();

        void calcShockFactor();

        /**\brief Time marching.*/
        virtual void marching() = 0;

        virtual void updateOld() {}

        /**\brief Update right hand side*/
        virtual void updateRhs() {}

        virtual void updateWStar() {}

        virtual hur_nodiscard bool isBDFScheme() const noexcept;

        virtual hur_nodiscard bool isESDIRKScheme() const noexcept;

        /** \brief Pseudo time step */
        pseudoTime &pseudoTimes();

        inline virtual hur_nodiscard bool updatedTemperature() const noexcept;
    };
} // namespace OpenHurricane

#include "timeMarching.inl"