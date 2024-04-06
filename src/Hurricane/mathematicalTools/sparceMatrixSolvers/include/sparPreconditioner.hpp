/*!
 * \file sparPreconditioner.hpp
 * \brief Header of Preconditioner.
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
#include "CRSMatrix.hpp"
#include "controller.hpp"
#include "objectFactory.hpp"

namespace OpenHurricane {
    namespace sparseSolver {

        template <class Type> class preconditioner {
        protected:
            const CRSMatrix<Type> &A_;

        public:
            using value_type = typename CRSMatrix<Type>::value_type;
            using solution_type = typename CRSMatrix<Type>::solution_type;

            declareClassName(preconditioner);

            declareObjFty(preconditioner, controller,
                          (const CRSMatrix<Type> &mat, const controller &cont), (mat, cont));

            inline preconditioner(const CRSMatrix<Type> &mat, const controller &cont) : A_(mat) {}

            hur_nodiscard static inline uniquePtr<preconditioner>
            creator(const CRSMatrix<Type> &mat, const controller &cont) {
                string preType = cont.findWord(preconditioner::className_);
                defineInObjCreatorTmpl(preconditioner, preType, controller, (mat, cont));
            }

            inline virtual ~preconditioner() noexcept {}

            /**
             * \brief Solve or smooth the solution for a given number of iterations.
             */
            virtual void precondition(Array<solution_type> &x,
                                      const Array<solution_type> &b) const = 0;
        };

        using realPreconditioner = preconditioner<real>;
        using realBlockPreconditioner = preconditioner<realSquareMatrix>;

    } // namespace sparseSolver
} // namespace OpenHurricane