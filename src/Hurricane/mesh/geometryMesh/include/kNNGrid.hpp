/*!
 * \file kNearestNeighbour.hpp
 * \brief Header of kNearestNeighbour.
 *       The subroutines and functions are in the <i>kNearestNeighbour.cpp</i> file.
 * \author Yang Hongzhen
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
#include "boundBox.hpp"

namespace OpenHurricane {

    /*!
     * \brief Base class of 3-Dimensional kNN algorithm.
     * Expansion to k-dimension is to be done.
     */
    class kNearestNeighbour {
    protected:
        /*!\brief Number of interpolation source base cells.*/
        integer k_;

        /*!\brief Source mesh cell centre vectorArray.*/
        const vectorArray &sor_;

        /*!\brief Target mesh cell centre vectorArray.*/
        const vectorArray &tar_;

        /*!\brief Nearest k elements of target cell.*/
        mutable uniquePtr<integerListList> nbrPtr_;

        virtual void getNearestNeigbhour() = 0;

    public:
        kNearestNeighbour() = delete;

        /*!\brief Construct from controller.*/
        kNearestNeighbour(const controller &cont, const vectorArray &sor, const vectorArray &tar);

        /*!\brief Construct from controller.*/
        kNearestNeighbour(const integer k, const vectorArray &sor, const vectorArray &tar);

        kNearestNeighbour(const kNearestNeighbour &) = delete;
        kNearestNeighbour &operator=(const kNearestNeighbour &) = delete;

        virtual ~kNearestNeighbour() noexcept { nbrPtr_.clear(); }

        inline integerListList &nearestNbr() {
            if (!nbrPtr_) {
                getNearestNeigbhour();
            }
            return *nbrPtr_;
        }

        hur_nodiscard inline integer k() const noexcept { return k_; }
    };

    /*!
     * \brief Derive class of 3-Dimensional kNN algorithm.
     * Based on grid, kdtree-based is to be done.
     */
    class kNNGrid : public kNearestNeighbour {
    private:
        /*!\brief BoundBox of cell centre.*/
        boundBox bb_;

        /*!\brief Length scale to divide bound box.*/
        vector l_;

        /*!\brief Length ref to normalize bound box and points.*/
        real ref_;

        /*!\brief bound box min point.*/
        List<integerListListList> contain_;

    protected:
        virtual void getNearestNeigbhour();

    private:
        hur_nodiscard integerVector getIndex(const vector point) const;

        void creatBoundBox();

        void expandIJK(integer &i1, integer &i2, integer &j1, integer &j2, integer &k1, integer &k2,
                       const integer nx, const integer ny, const integer nz, const realArray &uu,
                       const real dm, bool &flag, integerList &ib);

        void expandIJK(integerVector &ijk1, integerVector &ijk2, const integer nx, const integer ny,
                       const integer nz, const realArray &uu, const real dm, bool &flag,
                       integerList &ib);

    public:
        kNNGrid() = delete;

        /*!\brief Construct from controller.*/
        kNNGrid(const controller &cont, const vectorArray &sor, const vectorArray &tar);

        /*!\brief Construct from controller.*/
        kNNGrid(const integer k, const vectorArray &sor, const vectorArray &tar);

        kNNGrid(const kNearestNeighbour &) = delete;
        kNNGrid &operator=(const kNearestNeighbour &) = delete;

        inline virtual ~kNNGrid() noexcept {}
    };
} // namespace OpenHurricane