﻿/*!
 * \file TVDRK4.hpp
 * \brief Headers of class of the 4th order TVD Runge-Kutta.
 *        The subroutines and functions are in the <i>TVDRK4.cpp</i> file.
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

    /*!\brief The base class for 4th order TVD RK scheme.*/
    class TVDRK4 : public timeMarching {
    private:
        /**
         * \brief Coefs. at each stage.
         * i.e.
         *      1/2
         *      1/2
         *      1
         *      1/6
         */
        realArray alpha_;

        /**
         * \brief Coefs. for summing residuals in the last stage.
         * i.e.
         *      1
         *      2
         *      2
         *      1
         */
        realArray beta_;

        /**
         * \brief The number of stages.
         */
        integer stageSize_;

        /**
         * \brief The index of density in the list of primitive objects.
         */
        integer rhoOi_;

    public:
        declareClassNames;

        /*!\brief Disallow null constructor.*/
        TVDRK4() = delete;

        /*!\brief Disallow copy constructor.*/
        TVDRK4(const TVDRK4 &) = delete;

        /**
         * \brief Construct from components.
         */
        TVDRK4(const controller &cont, const runtimeMesh &mesh, const flowModel &flowM,
               solver &_solver, cellVectorArray &v);

        /*!\brief Destructor.*/
        virtual ~TVDRK4() noexcept;

        void stageUpdate();

        void updatePrimitives();

        void restoreOldResiduals(const integer stagei);
        void initOldResiduals();

        virtual void timeStep();

        virtual void initializing();

        virtual void marching();
    };

} // namespace OpenHurricane