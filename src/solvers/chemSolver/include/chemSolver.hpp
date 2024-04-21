/*!
 * \file chemSolver.hpp
 * \brief Headers of solver for chemical ODEs.
 *        The subroutines and functions are in the <i>chemSolver.cpp</i> file.
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

#include "chemistrySource.hpp"
#include "solver.hpp"

namespace OpenHurricane {
    /*!\brief The class for chemical ODEs solver.*/
    class chemSolver : public solver {
    protected:
        /*!\brief Delete the pointer.*/
        template <class Type> inline void deletePointer(geometryArray<Type, cellMesh> *_gPtr) const;

        uniquePtr<chemistrySource> chemistryPtr_;

        realArray yyi_;

        real p_;
        real rho_;
        real T_;

        real timeStep_;
        real endTime_;

        bool writeMolarFraction_;

        integer printStep_;
        integer writeStep_;

        integerList fuelSpecId_;

        /**
         * \brief The number of atoms C in species s.
         */
        integerList NCs_;

        /**
         * \brief The number of atoms H in species s.
         */
        integerList NHs_;

        /**
         * \brief The number of atoms O in species s.
         */
        integerList NOs_;

        /**
         * \brief The number of atoms N in species s.
         */
        integerList NNs_;

        /**
         * \brief The proportion offuel oxygen to fuel carbon.
         */
        real zp_;

        bool writePhiProgress_;

        void calcPhiAndPhiProgress(const real rho, const realArray &yi, real &phi, real &phip,
                                   real &phil) const;

        bool writeHeatReleaseRate_;

        bool writeSpeciesTimeScale_;

        real ignitTemThreshold_;

    public:
        declareClassName(chemSolver);

        chemSolver(iteration &iter, const runtimeMesh &mesh);

        /*!\brief Destructor.*/
        virtual ~chemSolver() noexcept;

        /**\brief Solve the problem.*/
        virtual void solve();

        virtual void solving();

        virtual void BDFSolve();

        /*!\brief Clear the solver.*/
        virtual void clear() noexcept;

        virtual void bc();

        virtual void updateProperties();

        virtual void initialize();

        void timeStep(realArray &dt) {}

        virtual void calculateFc();

        virtual void calculateFv() {}

        virtual void calculateSource() {}

        virtual void updatePrimitives(const bool shouldUpdateTemp = false);

        virtual void updateFlowOld() {}

        virtual void write();
    };
} // namespace OpenHurricane

#include "chemSolver.inl"