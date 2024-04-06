/*!
 * \file nozzleRegion.hpp
 * \brief Headers of the nozzle region of mesh.
 *        The subroutines and functions are in the <i>nozzleARegion.cpp</i> file.
 * \author Peng Jian
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
#include <map>
#include <string>

namespace OpenHurricane {
    class nozzleRegion : public markRegion {
    private:
        /*!\brief Axis max [Unit: m].*/
        vector Aexit_;

        /*!\brief Axis min [Unit: m].*/
        vector Ahead_;

        /*!\brief Axis throat [Unit: m].*/
        vector Athroat_;

        /*!\brief The diameter of the head [Unit: m].*/
        real Dhead_;

        /*!\brief The diameter of the exit [Unit: m].*/
        real Dexit_;

        /*!\brief The diameter of the throat [Unit: m].*/
        real Dthroat_;

        /*!\brief The specific heat ratio of rocket jet.*/
        real gammaConst_;

        /*!\brief The specific heat ratio of rocket jet.*/
        real criticalValue_;

        /*!\brief Combustion chamber pressure*/
        real pc_;

        /*!\brief Inlet boundary surface name.*/
        std::string boundaryFaceName_;

        /*!\brief The pressure has been patched in this nozzle region or not.*/
        bool pressPatchedFlag_;

        /*!\brief The flag map of varibles those should been specified.*/
        mutable uniquePtr<std::map<std::string, std::string>> varPatchFlagMap_;

        hur_nodiscard real divePresFromDia(const real diaRatio) const;

        hur_nodiscard real convPresFromDia(const real diaRatio) const;

        hur_nodiscard real areaRatioFromPres(const real presRatio) const;

        void calTemperature(realGeometryArray<cellMesh> &T) const;

        void calDensity(realGeometryArray<cellMesh> &rho) const;

        void calPressure(realGeometryArray<cellMesh> &p) const;

    public:
        declareClassNames;

        nozzleRegion(const controller &cont);

        /*!\brief Destructor.*/
        inline virtual ~nozzleRegion() noexcept { varPatchFlagMap_.clear(); }

        // Member functions

        hur_nodiscard inline const vector &Aexit() const noexcept;
        hur_nodiscard inline const vector &Ahead() const noexcept;
        hur_nodiscard inline const vector &Athroat() const noexcept;
        hur_nodiscard inline real Dhead() const noexcept;
        hur_nodiscard inline real Dexit() const noexcept;
        hur_nodiscard inline real Dthroat() const noexcept;
        hur_nodiscard inline real gammaConst() const noexcept;
        hur_nodiscard inline const std::string &bcfName() const noexcept;

        void setGammaConst(const real gammabc);

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

#include "nozzleRegion.inl"
