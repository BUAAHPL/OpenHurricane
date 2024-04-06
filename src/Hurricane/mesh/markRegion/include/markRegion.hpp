/*!
 * \file markRegion.hpp
 * \brief Headers of the mark region of mesh.
 *        The subroutines and functions are in the <i>markRegion.cpp</i> file.
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

#include "geometryMesh.hpp"
#include "geometryArrays.hpp"

namespace OpenHurricane {
    class cellMesh;
    class runtimeMesh;

    class markRegion {
    public:
        enum class regionOptions : short { INSIDE, OUTSIDE, NO_OPTIONs };

    private:
        /*!\brief Index of the region.*/
        integer id_;

        /*!\brief The option of the region.*/
        regionOptions option_;

    public:
        declareClassNames;
        declareObjFty(markRegion, controller, (const controller &cont), (cont));

        inline markRegion();

        /*!\brief Construct form controller.*/
        markRegion(const controller &cont);

        static uniquePtr<markRegion> creator(const controller &cont);

        /*!\brief Destructor.*/
        virtual ~markRegion() noexcept {}

        /*!\brief Return the index of this region.*/
        hur_nodiscard inline integer id() const noexcept;

        /*!\brief Return true if the option set.*/
        hur_nodiscard inline bool isOptionSet() const noexcept;

        /*!\brief Return true if the option set.*/
        hur_nodiscard inline regionOptions option() const noexcept;

        /*!\brief Return true if the option is inside.*/
        hur_nodiscard inline bool isInside() const noexcept;

        /*!\brief Return true if the option is outside.*/
        hur_nodiscard inline bool isOutside() const noexcept;

        // interface

        hur_nodiscard virtual integerList regionCellId(const runtimeMesh &mesh) const;

        /*!
         * \brief To patch different values of flow variables into different cells.
         * \return True if patching success.
         */
        virtual bool patching(realGeometryArray<cellMesh> &cellQ, real &value) const = 0;

        /*!
         * \brief To patch different values of flow variables into different cells.
         * \return True if patching success.
         */
        virtual bool patching(vectorGeometryArray<cellMesh> &cellQ, vector &value) const = 0;

        
        /*!
         * \brief To patch different values of flow variables into different cells.
         * \return True if patching success.
         */
        virtual bool distributing(realGeometryArray<cellMesh> &cellQ, std::string &value) const = 0;

        /*!
         * \brief To patch different values of flow variables into different cells.
         * \return True if patching success.
         */
        virtual bool distributing(vectorGeometryArray<cellMesh> &cellQ,
                                  std::string &value) const = 0;
    };

} // namespace OpenHurricane

#include "markRegion.inl"