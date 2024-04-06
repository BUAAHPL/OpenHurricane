/*!
 * \file thirdUnsMUSCL.hpp
 * \brief Header of third order unstructured MUSCL interpolation.
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
    class thirdUnsMUSCL : public reconstruction {
    private:
        inline real beta(const real r) const;
        vector beta(const vector &r) const;

        inline real phi(const real r) const;
        vector phi(const vector &r) const;

    public:
        declareClassNames;

        /*!\brief Construct as null.*/
        thirdUnsMUSCL();

        /*!\brief Construct from controller.*/
        thirdUnsMUSCL(const controller &cont);

        /*!\brief Disallow copy constructor.*/
        thirdUnsMUSCL(const thirdUnsMUSCL &) = delete;

        /*!\brief Disallow bitwise assignment.*/
        thirdUnsMUSCL &operator=(const thirdUnsMUSCL &) = delete;

        /*!\brief Destructor.*/
        virtual ~thirdUnsMUSCL() noexcept;

        /**
         * \brief  MUSCL reconstruction
         * \param[in] QCL The face left side cell variable
         * \param[in] QCR The face right side cell variable
         * \param[in] gradCL The face left side cell gradient
         * \param[in] gradCR The face right side cell gradient
         * \param[in] dLR The distance vector from left cell to right cell
         *                Use with care: dLR = cellCentre(r) - cellCentre(l)
         * \param[out] ql The face left side variable
         * \param[out] qr The face right side variable
         */
        void UnsMUSCLRec(const real QCL, const real QCR, const vector &gradCL, const vector &gradCR,
                         const vector &dLR, real &ql, real &qr) const;

        /**
         * \brief  MUSCL reconstruction
         * \param[in] QCL The face left side cell variable
         * \param[in] QCR The face right side cell variable
         * \param[in] gradCL The face left side cell gradient
         * \param[in] gradCR The face right side cell gradient
         * \param[in] dLR The distance vector from left cell to right cell
         *                Use with care: dLR = cellCentre(r) - cellCentre(l)
         * \param[out] ql The face left side variable
         * \param[out] qr The face right side variable
         */
        void UnsMUSCLRec(const vector QCL, const vector QCR,
                         const typename outerProduct<vector, vector>::type &gradCL,
                         const typename outerProduct<vector, vector>::type &gradCR,
                         const vector &dLR, vector &ql, vector &qr) const;

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

#include "thirdUnsMUSCL.inl"