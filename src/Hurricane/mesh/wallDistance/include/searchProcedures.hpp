/*!
 * \file searchProcedures.hpp
 * \brief Headers of class of search procedures of distance method.
 *        The subroutines and functions are in the <i>searchProcedures.cpp</i> file.
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

#include "distanceMethod.hpp"

namespace OpenHurricane {

    /*!\brief face zone package: face centre and face area.*/
    class faceZonePackage {
    private:
        /*!\brief Face centre.*/
        vectorArray faceCentre_;

        /*!\brief Face area.*/
        vectorArray faceArea_;

        vectorArrayArray fp_;

    public:
        faceZonePackage();

        /*!\brief Destructor.*/
        ~faceZonePackage() noexcept {}

        /*!\brief Face centre.*/
        hur_nodiscard inline vectorArray &faceCentre() noexcept;

        /*!\brief Face centre.*/
        hur_nodiscard inline const vectorArray &faceCentre() const noexcept;

        /*!\brief Face area.*/
        hur_nodiscard inline vectorArray &faceArea() noexcept;

        /*!\brief Face area.*/
        hur_nodiscard inline const vectorArray &faceArea() const noexcept;

        /*!\brief Face point.*/
        hur_nodiscard inline vectorArrayArray &fp() noexcept;

        /*!\brief Face point.*/
        hur_nodiscard inline const vectorArrayArray &fp() const noexcept;

        void clear() noexcept;
    };

    /*!\brief The class of search procedure for distance calculation.*/
    class searchProcedures : public distanceMethod {
    private:
        /*!\brief Zone package.*/
        List<faceZonePackage> facePack_;

        integerList faceIdList_;

        integerList iw1_;

        integerList iw2_;

        integerVector2DList link_;

        integerList itt_;

        integer nFaces_;

        /*!\brief Should normal to the wall distance.*/
        bool normal_;

        /*!\brief Broadcast face packages.*/
        void setFacePack();

        void setFacePoint();

        void distTransfer();
        void smooth(cellRealArray &y, integer &t);

        real calcDist(const integer ip, const integer ii, const vector &cC) const;

    public:
        declareClassNames;

        searchProcedures() = delete;

        searchProcedures(const controller &cont, const runtimeMesh &mesh,
                         const integerList zoneIdL);

        searchProcedures(const runtimeMesh &mesh, const integerList zoneIdL);

        searchProcedures(const searchProcedures &) = delete;
        searchProcedures &operator=(const searchProcedures &) = delete;

        inline virtual ~searchProcedures() noexcept {}

        virtual bool getDistance(cellRealArray &y);

        virtual bool getDistance(cellRealArray &y, cellVectorArray &ny);
    };
} // namespace OpenHurricane

#include "searchProcedures.inl"