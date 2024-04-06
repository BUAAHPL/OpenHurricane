/*!
 * \file boundary.hpp
 * \brief Headers of the boundary.
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
#include "geometryArray.hpp"
#include "meshElements.hpp"
#include "objectFactory.hpp"
#include "registerTable.hpp"
#include "runtimeMesh.hpp"
#include "smartPointerList.hpp"
#include "zoneMesh.hpp"

namespace OpenHurricane {

    template <class Type, class geoMesh> class boundary {
    protected:
        /*!\brief Reference to faceZone.*/
        const faceZone &boundaryZone_;

        /*!\brief Reference to geometryArray.*/
        geometryArray<Type, geoMesh> &varArray_;

        bool specified_;

        /*!\brief Initialization value with unit.*/
        mutable uniquePtr<Type> initValue_;

    public:
        declareClassNames;
        declareObjFty(boundary, controller,
                      (const faceZone &fZ, geometryArray<Type, geoMesh> &gM,
                       const controller &cont),
                      (fZ, gM, cont));

        boundary(const faceZone &fZ, geometryArray<Type, geoMesh> &gf, const controller &cont)
            : boundaryZone_(fZ), varArray_(gf), specified_(false), initValue_(nullptr) {

            const controller &initCont = cont.topCont().subController("initialization");
            if (initCont.findWord("initFromBoundary") == boundaryZone_.name()) {
                initValue_.reset(new Type());
                if (cont.found(varArray_.name())) {
                    readValue(*initValue_, varArray_.name(), cont);
                } else {
                    *initValue_ = Zero;
                }
            }
        }

        /**\brief Construct from components.*/
        boundary(const boundary<Type, geoMesh> &bB)
            : boundaryZone_(bB.boundaryZone_), varArray_(bB.varArray_), specified_(bB.specified_),
              initValue_(nullptr) {
            if (bB.initValue_) {
                initValue_.reset(new Type());
                *initValue_ = *bB.initValue_;
            }
        }

        /*!\brief Return a pointer to a new baseBoundary created on freestore given
        faceZone and geometryArray.*/
        hur_nodiscard static uniquePtr<boundary<Type, geoMesh>>
        creator(const faceZone &fZ, geometryArray<Type, geoMesh> &gf, const controller &cont) {
            const controller &fZCont =
                cont.subController("boundaryCondition").subController(fZ.name());
            string bcType;
            if (fZCont.found("defultType")) {
                if (fZCont.found(gf.name() + string("bcType"))) {
                    bcType = fZCont.findWord(gf.name() + string("bcType"));
                    typename controllercreactorMapType::iterator ptrIter =
                        controllercreatorMapPtr_->find(static_cast<std::string>(bcType));
                    if (ptrIter == controllercreatorMapPtr_->end()) {
                        LFatal("Unknown boundary condiction: %s \nValid boundary condiction of "
                               "current "
                               "program are: ",
                               bcType.c_str(), stringMapDoc(*controllercreatorMapPtr_).c_str());
                    }
                    return (ptrIter->second)(fZ, gf, fZCont);
                } else {
                    bcType = fZCont.findWord("defultType");
                }
            } else {
                bcType = fZCont.findWord("bcType");
            }
            typename controllercreactorMapType::iterator ptrIter =
                controllercreatorMapPtr_->find(static_cast<std::string>(bcType));

            if (ptrIter == controllercreatorMapPtr_->end()) {
                LFatal("Unknown boundary condiction: %s \nValid boundary condiction of current "
                       "program are: ",
                       bcType.c_str(), stringMapDoc(*controllercreatorMapPtr_).c_str());
            }

            return (ptrIter->second)(fZ, gf, fZCont);
        }

        virtual ~boundary() noexcept {}

        /*!\brief Return boundary faceZone.*/
        hur_nodiscard inline const faceZone &boundaryfZ() const { return boundaryZone_; }

        hur_nodiscard inline const geometryArray<Type, geoMesh> &varField() const noexcept {
            return varArray_;
        }

        hur_nodiscard inline geometryArray<Type, geoMesh> &varField() noexcept { return varArray_; }

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept { return varArray_.mesh(); }

        /*!\breif Return initValue.*/
        hur_nodiscard inline const Type &initValue() const {
            if (!initValue_) {
                initValue_.reset(new Type());
                *initValue_ = Zero;
            }
            return *initValue_;
        }

        virtual void updateBoundary() = 0;

        void setSpecified() noexcept { specified_ = true; }
        void readValue(Type &value, const std::string &name, const controller &cont) const {
            LFatal("This type has not been specified in readValue(...) function.");
        }

        void readVector(vector &value, const std::string &name, const controller &cont) const {
            LFatal("This type has not been specified in readValue(...) function.");
        }
    };
} // namespace OpenHurricane