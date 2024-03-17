/*!
 * \file packs.hpp
 * \brief Header of packs
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
#include "ArrayInclude.hpp"
#include "HurMPIBase.hpp"
#include "object.hpp"
#include "Lists.hpp"

namespace OpenHurricane {
    template <class Type> using packs = Array<Type>;

    using scalarPacks = packs<real>;

    using blockPacks = packs<Array<real>>;

    class packageMap {
    protected:
        List<object *> objectList_;

        /**
         * \brief The number of packgae scalar parameters.
         */
        integer countPackage_;

        integerVector2DList packageMap_;

    public:
        inline packageMap() : objectList_(), countPackage_(0), packageMap_() {}

        inline virtual ~packageMap() noexcept { clear(); }

        void unbind() noexcept;
        void clear() noexcept;

        void addObject(object &ob);

        /**
         * \brief The number of packgae scalar parameters.
         */
        hur_nodiscard inline integer countPackage() const noexcept { return countPackage_; }
    };

} // namespace OpenHurricane
