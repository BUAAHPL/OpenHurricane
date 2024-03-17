/*!
 * \file fixedGradient.hpp
 * \brief Headers of the fixedGradient.
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
#include "boundaries.hpp"

namespace OpenHurricane {

    template <class Type> class fixedGradient : public boundary<Type, cellMesh> {
    public:
        using Base = boundary<Type, cellMesh>;

    private:
        /*!\brief Fixed normal gradient of boundary.*/
        Type gradient_;

    public:
        declareClassNames;

        fixedGradient(const faceZone &fZ, geometryArray<Type, cellMesh> &gf,
                       const controller &cont)
            : Base(fZ, gf, cont) {
            this->readValue(gradient_, gf.name(), cont);
            this->setSpecified();
        }

        fixedGradient(const fixedGradient &bB) : Base(bB), gradient_(bB.gradient_) {}

        virtual ~fixedGradient() noexcept {}

        virtual void updateBoundary() {
            auto &var = Base::varArray_;
            const faceList &fL = var.mesh().faces();
            const auto &fa = var.mesh().faceArea();
            const auto &cc = var.mesh().cellCentre();

            for (integer fi = this->boundaryZone_.firstIndex();
                 fi < this->boundaryZone_.lastIndex() + 1; fi++) {
                const auto cl = fL[fi].leftCell();
                const auto cr = fL[fi].rightCell();
                real delta = (cc[cr] - cc[cl]) * fa[fi].normalized();
                var[cr] = var[cl] + delta * gradient_;
            }
        }
    };

    using realFixedGradient = fixedGradient<real>;
    using vectorFixedGradient = fixedGradient<vector>;
} // namespace OpenHurricane