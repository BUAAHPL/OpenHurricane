/*!
 * \file Linear.hpp
 * \brief Header of piecewise linear reconstruction.
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

#include "limitersForLinear.hpp"
#include "reconstruction.hpp"

namespace OpenHurricane {

    class Linear : public reconstruction {
    private:
        /*!\brief The pointer to the limiter function.*/
        uniquePtr<limitersForLinear> limiterPtr_;

    public:
        declareClassNames;

        /*!\brief Construct as null.*/
        Linear();

        /*!\brief Construct from controller.*/
        Linear(const controller &cont);

        /*!\brief Destructor.*/
        virtual ~Linear() noexcept;

        inline virtual void calcLimiter(geometryArray<real, cellMesh> &cellQ) const;
        inline virtual void calcLimiter(geometryArray<vector, cellMesh> &cellQ) const;

        /**
         * \brief piecewise linear reconstruction
         * \param[in] QCL The face left side cell variable
         * \param[in] QCR The face right side cell variable
         * \param[in] gradCL The face left side cell gradient
         * \param[in] gradCR The face right side cell gradient
         * \param[in] rL The distance vector from left cell centroid to face-midpoint
         * \param[in] rR The distance vector from right cell centroid to face-midpoint
         * \param[out] ql The face left side variable
         * \param[out] qr The face right side variable
         */
        template <class Type>
        void linearRecon(const Type &QCL, const Type &QCR,
                         const typename outerProduct<vector, Type>::type &gradCL,
                         const typename outerProduct<vector, Type>::type &gradCR, const vector &rL,
                         const vector &rR, const Type &limiterL, const Type &limiterR, Type &ql,
                         Type &qr) const;

        virtual void calcReconstruction(const geometryArray<real, cellMesh> &cellQ,
                                        const integer faceI, real &ql, real &qr) const;

        virtual void calcReconstruction(const geometryArray<vector, cellMesh> &cellQ,
                                        const integer faceI, vector &ql, vector &qr) const;

        virtual void calcReconstruction(const geometryArray<real, cellMesh> &cellQ,
                                        const integer faceZoneI, realArray &ql,
                                        realArray &qr) const;

        virtual void calcReconstruction(const geometryArray<vector, cellMesh> &cellQ,
                                        const integer faceZoneI, vectorArray &ql,
                                        vectorArray &qr) const;

        virtual void calcReconstructionWithoutLimit(const geometryArray<real, cellMesh> &cellQ,
                                                    const integer faceI, real &ql, real &qr) const;

        virtual void calcReconstructionWithoutLimit(const geometryArray<vector, cellMesh> &cellQ,
                                                    const integer faceI, vector &ql,
                                                    vector &qr) const;

        virtual void calcReconstructionWithoutLimit(const geometryArray<real, cellMesh> &cellQ,
                                                    const integer faceZoneI, realArray &ql,
                                                    realArray &qr) const;

        virtual void calcReconstructionWithoutLimit(const geometryArray<vector, cellMesh> &cellQ,
                                                    const integer faceZoneI, vectorArray &ql,
                                                    vectorArray &qr) const;
    };

} // namespace OpenHurricane

#include "Linear.inl"