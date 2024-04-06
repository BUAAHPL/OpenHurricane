/*!
 * \file RBF.hpp
 * \brief Header of Radial Basis Function (RBF).
 *       The subroutines and functions are in the <i>RBF.cpp</i> file.
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

#include "OpenHurricane.hpp"
#include "compactSupportFunction.hpp"
#include "matrices.hpp"
#include "meshElements.hpp"
#include "zoneMesh.hpp"

namespace OpenHurricane {
    class geometryMesh;

    /**
     * \brief The base class of RBF.
     */
    class RBF {
    private:
        const geometryMesh &mesh_;

        const faceZone &movingZone_;

        /*!\brief Basic points for mesh motion. It is same in all processes.*/
        pointField basePoints_;

        pointField oldBasePoint_;

        uniquePtr<compactSupportFunction> CSFPtr_;

        /*!\brief Distance matrix.*/
        realSymmMatrix M_;

        /*!\brief Interpolation coefficients for coordinate x.*/
        realArray alphax_;

        /*!\brief Interpolation coefficients for coordinate y.*/
        realArray alphay_;

        /*!\brief Interpolation coefficients for coordinate z.*/
        realArray alphaz_;

        /*!\brief Selected boundary points label list.*/
        integerList selectedPointList_;

        bool firstGuess_;

        /*!\brief Tolerance.*/
        real epsilon_;

    public:
        declareClassName(RBF);

        RBF() = delete;
        RBF(const RBF &) = delete;
        RBF &operator=(const RBF &) = delete;

        /*!\brief Construct from controller.*/
        RBF(const controller &cont, const geometryMesh &mesh, const faceZone &movingZone,
            const integer fzid);

        inline virtual ~RBF() noexcept { CSFPtr_.clear(); }

        void setBasePoints(const vectorArray &bp);
        void setOldBasePoints();

        hur_nodiscard inline real Derror(const vector &x0, const vector &x1) const;

        /*!\brief Offset vector for point pi.*/
        hur_nodiscard vector g(const vector &pi) const;

        /*!\brief Offset vector for point pi.*/
        void g(const vector &pi, vector &gg) const;

        /*!\brief New point for point pi.*/
        hur_nodiscard vector newPoint(const vector &oldPi) const;

        /*!\brief New point for point pi.*/
        void newPoint(const vector &oldPi, vector &np) const;

        void getAlpha();

        void initialSelected();

        void solving();

        hur_nodiscard bool isSatisfied(const real err) const;
    };

} // namespace OpenHurricane

#include "RBF.inl"