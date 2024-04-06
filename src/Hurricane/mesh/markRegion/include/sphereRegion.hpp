/*!
 * \file sphereRegion.hpp
 * \brief Headers of the sphere region of mesh.
 *        The subroutines and functions are in the <i>sphereRegion.cpp</i> file.
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
#include "markRegion.hpp"

namespace OpenHurricane {
    class sphereRegion : public markRegion {
    private:
        /*!\brief The center of sphere [Unit: m]*/
        vector center_;

        /*!\brief The radius of the sphere region [Unit: m].*/
        real radius_;

    public:
        declareClassNames;

        sphereRegion(const controller &cont);

        /*!\brief Destructor.*/
        inline virtual ~sphereRegion() noexcept {}

        hur_nodiscard inline vector center() const noexcept;
        hur_nodiscard inline real radius() const noexcept;

        hur_nodiscard virtual integerList regionCellId(const runtimeMesh &mesh) const;

        /*!
         * \brief To patch different values of flow variables into different cells.
         * \param[in] cellQ - The field variable.
         * \param[in] value - The value to be patched for the given variable. Must given with units.
         * \return True if patching success.
         */
        virtual bool patching(realGeometryArray<cellMesh> &cellQ, real &value) const;

        /*!
         * \brief To patch different values of flow variables into different cells.
         * \param[in] cellQ - The field variable.
         * \param[in] value - The value to be patched for the given variable. Must given with units.
         * \return True if patching success.
         */
        virtual bool patching(vectorGeometryArray<cellMesh> &cellQ, vector &value) const;

        /*!
         * \brief To patch different values of flow variables into different cells.
         * \param[in] cellQ - The field variable.
         * \param[in] value - The value to be distributed for the given formula expression. Must given with units.
         * \return True if distributing success.
         */
        virtual bool distributing(realGeometryArray<cellMesh> &cellQ, std::string &value) const;

        /*!
         * \brief To patch different values of flow variables into different cells.
         * \param[in] cellQ - The field variable.
         * \param[in] value - The value to be distributed for the given formula expression. Must given with units.
         * \return True if distributing success.
         */
        virtual bool distributing(vectorGeometryArray<cellMesh> &cellQ, std::string &value) const;
    private:
        /*!
         * \brief To patch different values of flow variables into different cells.
         * \param[in] cellQ - The field variable.
         * \param[in] value - The value to be patched for the given variable. Must given with units.
         * \return True if patching success.
         */
        template <class Type, class meshType>
        bool patch(geometryArray<Type, meshType> &cellQ, Type &value) const;
    };
} // namespace OpenHurricane

#include "sphereRegion.inl"