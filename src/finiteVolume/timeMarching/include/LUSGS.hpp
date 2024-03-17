/*!
 * \file LUSGS.hpp
 * \brief Headers of class of the LUSGS.
 *        The subroutines and functions are in the <i>LUSGS.cpp</i> file.
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
#include "timeMarching.hpp"

namespace OpenHurricane {

    /*!\brief The base class for LUSGS scheme.*/
    class LUSGS : public timeMarching {
    private:
        real omega_;

        real betaT_;

    protected:
        bool modifiedDiagnalTime_;

        real alphaDT_;

        real scaleTimeStep_;

        boolList *fullJacFlagPtr_;

        hur_nodiscard inline boolList &fullJacFlag();

    public:
        declareClassNames;

        /*!\brief Disallow null constructor.*/
        LUSGS() = delete;

        /*!\brief Disallow copy constructor.*/
        LUSGS(const LUSGS &) = delete;

        LUSGS(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
              solver &_solver, cellVectorArray &v, const bool modifiedDiagnalTime = false);

        /*!\brief Destructor.*/
        virtual ~LUSGS() noexcept;

        void lowerLoop();
        void upperLoop();
        void update();

        void fluxJacoDQ(const integer i, const vector &n, realArray &adq);

        virtual void marching();

    protected:
    public:
        static void lowerLoop(cellRealArray &_cellQ, const vector2DArray &_ra, const realArray &_dt,
                              realArray &_dq);
        static void upperLoop(cellRealArray &_cellQ, const vector2DArray &_ra, const realArray &_dt,
                              realArray &_dq);

        inline virtual hur_nodiscard bool updatedTemperature() const noexcept;
    };

} // namespace OpenHurricane

#include "LUSGS.inl"