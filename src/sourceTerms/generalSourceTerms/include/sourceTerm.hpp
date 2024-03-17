/*!
 * \file sourceTerm.hpp
 * \brief Headers of base class of source term.
 *        The subroutines and functions are in the <i>sourceTerm.cpp</i> file.
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

#include "flowModel.hpp"
#include "geometryArray.hpp"
#include "iteration.hpp"

namespace OpenHurricane {

    /*!\brief The class of sourceTerm.*/
    class sourceTerm {
    protected:
        /**\brief The flow model.*/
        const flowModel &flows_;

        const iteration &iter_;

    public:
declareClassNames;
declareObjFty(sourceTerm,controller,
                            (const flowModel &_flows, const iteration &_iter,
                             const controller &_cont),
                            (_flows, _iter, _cont));

        // Constructors

        sourceTerm() = delete;

        sourceTerm(const flowModel &flows, const iteration &iter, const controller &cont);

        sourceTerm(const sourceTerm &st) = delete;

        static uniquePtr<sourceTerm> creator(const flowModel &flows, const iteration &iter,
                                           const controller &cont);

        static uniquePtr<sourceTerm> creator(const flowModel &flows, const iteration &iter,
                                           const controller &cont, const string &sourceType);

        /*!\brief Destructor.*/
        virtual ~sourceTerm() noexcept {}

        /**\brief The flow model.*/
        hur_nodiscard inline const flowModel &flows() const noexcept;

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept;
        hur_nodiscard inline const iteration &iter() const noexcept;

        /**
         * \brief Add to Continuity equation.
         */
        virtual inline void addSourceTerms(cellRealArray &rho) const;

        /**
         * \brief Add to momentum equation.
         */
        virtual inline void addSourceTerms(const cellRealArray &rho, cellVectorArray &U) const;

        /**
         * \brief Add to scalar equation.
         */
        virtual inline void addSourceTerms(const cellRealArray &rho, cellRealArray &phi) const;

        void operator=(const sourceTerm &st) = delete;
    };
} // namespace OpenHurricane

#include "sourceTerm.inl"