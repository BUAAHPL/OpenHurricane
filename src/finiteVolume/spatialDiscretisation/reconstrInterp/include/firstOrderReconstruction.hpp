/*!
 * \file firstOrderReconstruction.hpp
 * \brief Header of first order reconstruction.
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

#include "reconstruction.hpp"

namespace OpenHurricane {
    class firstOrderReconstruction : public reconstruction {
    public:
        declareClassNames;

        /*!\brief Construct as null.*/
        firstOrderReconstruction();

        /*!\brief Construct from controller.*/
        firstOrderReconstruction(const controller &cont);

        /*!\brief Destructor.*/
        virtual ~firstOrderReconstruction() noexcept;

        inline void firstOrderRec(const real QCL, const real QCR, real &ql, real &qr) const;

        inline void firstOrderRec(const vector QCL, const vector QCR, vector &ql, vector &qr) const;

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
    };
} // namespace OpenHurricane

#include "firstOrderReconstruction.inl"