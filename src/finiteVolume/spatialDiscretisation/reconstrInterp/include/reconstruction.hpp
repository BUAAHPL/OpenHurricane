/*!
 * \file reconstruction.hpp
 * \brief Headers of the reconstruction.
 *        The subroutines and functions are in the <i>reconstruction.cpp</i> file.
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
#include "cellArrays.hpp"
#include "commonInclude.hpp"
#include "faceArrays.hpp"
#include "gradient.hpp"
#include "smartPointerList.hpp"

namespace OpenHurricane {
    /*!\brief The base class of reconstruction.*/
    class reconstruction {
    private:
        /*!\brief The pointer for gradient method.*/
        mutable uniquePtr<gradient> gradPtr_;

    public:
        declareClassName(reconstruction);
        declareObjFty(reconstruction, controller, (const controller &cont), (cont));

        /*!\brief Construct as null.*/
        reconstruction();

        /*!\brief Construct from controller.*/
        reconstruction(const controller &cont);

        static uniquePtr<reconstruction> creator(const controller &cont);

        /*!\brief Destructor.*/
        virtual ~reconstruction() noexcept;

        /*!\brief Return the gradient method.*/
        gradient &grad() const;

        inline virtual void calcGrad(geometryArray<real, cellMesh> &cellQ) const;

        inline virtual void calcGrad(geometryArray<vector, cellMesh> &cellQ) const;

        inline void calcGrad(PtrList<geometryArray<real, cellMesh>> &cellQList) const;
        inline void calcGrad(PtrList<geometryArray<vector, cellMesh>> &cellQList) const;

        inline virtual void calcLimiter(geometryArray<real, cellMesh> &cellQ) const;
        inline virtual void calcLimiter(geometryArray<vector, cellMesh> &cellQ) const;

        template <class Type>
        inline void calcLimiter(PtrList<geometryArray<Type, cellMesh>> &cellQ) const;

        virtual void calcReconstruction(const geometryArray<real, cellMesh> &cellQ,
                                        const integer faceI, real &ql, real &qr) const = 0;

        virtual void calcReconstruction(const geometryArray<vector, cellMesh> &cellQ,
                                        const integer faceI, vector &ql, vector &qr) const = 0;

        virtual void calcReconstruction(const geometryArray<real, cellMesh> &cellQ,
                                        const integer faceZoneI, realArray &ql,
                                        realArray &qr) const = 0;

        virtual void calcReconstruction(const geometryArray<vector, cellMesh> &cellQ,
                                        const integer faceZoneI, vectorArray &ql,
                                        vectorArray &qr) const = 0;

        virtual void calcReconstructionWithoutLimit(const geometryArray<real, cellMesh> &cellQ,
                                                    const integer faceI, real &ql, real &qr) const {
            calcReconstruction(cellQ, faceI, ql, qr);
        }

        virtual void calcReconstructionWithoutLimit(const geometryArray<vector, cellMesh> &cellQ,
                                                    const integer faceI, vector &ql,
                                                    vector &qr) const {
            calcReconstruction(cellQ, faceI, ql, qr);
        }

        virtual void calcReconstructionWithoutLimit(const geometryArray<real, cellMesh> &cellQ,
                                                    const integer faceZoneI, realArray &ql,
                                                    realArray &qr) const {
            calcReconstruction(cellQ, faceZoneI, ql, qr);
        }

        virtual void calcReconstructionWithoutLimit(const geometryArray<vector, cellMesh> &cellQ,
                                                    const integer faceZoneI, vectorArray &ql,
                                                    vectorArray &qr) const {
            calcReconstruction(cellQ, faceZoneI, ql, qr);
        }
    };
} // namespace OpenHurricane

#include "reconstruction.inl"