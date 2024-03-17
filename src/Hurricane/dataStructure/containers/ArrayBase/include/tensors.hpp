/*!
 * \file tensors.hpp
 * \brief Header of tensors
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

#include "basicFunctions.hpp"
#include "diagTensorTmpl.hpp"
#include "sphericalTensorTmpl.hpp"
#include "string.hpp"
#include "symmTensorTmpl.hpp"
#include "tensorTmpl.hpp"
#include "vectors.hpp"

namespace OpenHurricane {

    using tensor = Tensor<real>;
    using integerTensor = Tensor<integer>;
    using floatTensor = Tensor<floatReal>;
    using symmTensor = SymmTensor<real>;
    using integerSymmTensor = SymmTensor<integer>;
    using diagTensor = DiagTensor<real>;
    using integerDiagTensor = DiagTensor<integer>;
    using sphericalTensor = SphericalTensor<real>;
    using integerSphericalTensor = SphericalTensor<integer>;

    // Global Identity tensor
    static const Identity<real> I;

    template <> hur_nodiscard std::string toString(const tensor &st);
    template <> hur_nodiscard std::string toString(const integerTensor &st);

#ifdef HURRICANE_DP
    template <> hur_nodiscard std::string toString(const floatTensor &st);
#endif //HURRICANE_DP
    template <> hur_nodiscard std::string toString(const symmTensor &st);
    template <> hur_nodiscard std::string toString(const integerSymmTensor &st);
    template <> hur_nodiscard std::string toString(const diagTensor &st);
    template <> hur_nodiscard std::string toString(const integerDiagTensor &st);
    template <> hur_nodiscard std::string toString(const sphericalTensor &st);
    template <> hur_nodiscard std::string toString(const integerSphericalTensor &st);

    /**
     * \brief Kronecker symbol.
     * \param[in] i - 0~2
     * \param[in] j - 0~2
     */
    hur_nodiscard real deltaij(const integer i, const integer j);

    /**
     * \brief Levi-Cevita symbol.
     * \param[in] i - 0~2
     * \param[in] j - 0~2
     * \param[in] k - 0~2
     */
    hur_nodiscard real epsilonijk(const integer i, const integer j, const integer k);
} // namespace OpenHurricane
