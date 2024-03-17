/*!
 * \file CFL.hpp
 * \brief Header of CFL (Courant�CFriedrichs�CLevy) number.
 *       The subroutines and functions are in the <i>CFL.cpp</i> file.
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

#include "dataStructure.hpp"
#include "geometryArrays.hpp"
#include "objectFactory.hpp"

namespace OpenHurricane {

    class iteration;
    class runtimeMesh;

    /*!\brief The base class of CFL (Courant�CFriedrichs�CLevy) number.*/
    class CFL {
    protected:
        const iteration &iter_;
        const runtimeMesh &mesh_;

        /*!\brief The maximum CFL number.*/
        real CFLMax_;

        /*!\brief The minmum CFL number.*/
        real CFLMin_;

        /*!\brief The current used CFL number.*/
        mutable real cCFL_;

        mutable bool notChange_;

        /*!\brief The step interval for printing CFL.*/
        integer stepForPrintCFL_;

        /** \brief The steps that stand for the initial stage of caculation. */
        integer initialStage_;

    public:
        declareClassNames;

        declareObjFty(CFL, controller,
                      (const iteration &iter, const runtimeMesh &mesh, const controller &cont),
                      (iter, mesh, cont));

        /*!\brief Disallow null constructor.*/
        CFL() = delete;

        /*!\brief Construct from controller.*/
        CFL(const iteration &iter, const runtimeMesh &mesh, const controller &cont);

        /*!\brief Disallow copy constructor.*/
        CFL(const CFL &cfl) = delete;

        hur_nodiscard static uniquePtr<CFL> creator(const iteration &iter, const runtimeMesh &mesh,
                                                    const controller &cont);

        /*!\brief Destructor.*/
        virtual ~CFL() noexcept {}

        hur_nodiscard const iteration &iter() const noexcept;

        hur_nodiscard const runtimeMesh &mesh() const noexcept;

        /*!\brief The maximum CFL number.*/
        hur_nodiscard inline real CFLMax() const noexcept { return CFLMax_; }

        /*!\brief The maximum CFL number.*/
        hur_nodiscard inline real &CFLMax() noexcept { return CFLMax_; }

        /*!\brief The minmum CFL number.*/
        hur_nodiscard inline real CFLMin() const noexcept { return CFLMin_; }

        /*!\brief The minmum CFL number.*/
        hur_nodiscard inline real &CFLMin() noexcept { return CFLMin_; }

        inline void setCFLMax(const real newCFLMax) noexcept {
            CFLMax_ = newCFLMax;
            CFLMin_ = min(CFLMin_, newCFLMax);
            cCFL_ = min(cCFL_, CFLMax_);
            cCFL_ = max(cCFL_, CFLMin_);
        }

        inline void setCFLMin(const real newCFLMin) noexcept {
            CFLMin_ = newCFLMin;
            CFLMax_ = max(CFLMax_, newCFLMin);
            cCFL_ = min(cCFL_, CFLMax_);
            cCFL_ = max(cCFL_, CFLMin_);
        }

        virtual integer CFLConst() const noexcept;
        virtual integer CFLLinear() const noexcept;

        /*!\brief Return the CFL number for current step.
         * \return The CFL number for current step.
         */
        virtual real getCFL() const = 0;

        virtual void setRelativeCorrection(const real relCorre) const;

        virtual void setBreakdowns() const;

        void setNotChange() const;

        /*!\brief Print CFL number */
        void printCFL() const;

        hur_nodiscard bool isInitialStage() const noexcept;
    };

} // namespace OpenHurricane