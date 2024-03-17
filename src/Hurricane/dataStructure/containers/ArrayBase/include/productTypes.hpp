/*!
 * \file vectorSpace.hpp
 * \brief Header of vector
 *       The subroutines and functions are in the <i>vectorSpace.cpp</i> file.
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
#include "feature.hpp"

namespace OpenHurricane {

    /*!\brief The abstract template class of the inner-product of two VectorSpace.*/
    template <class Type, class VectorSpace1, class VectorSpac2> class innerProductType {};

    /*!\brief The abstract template class of the outer-product of two VectorSpace.*/
    template <class Type, class VectorSpace1, class VectorSpace2> class outerProductType {};

    /*!\brief The abstract template class of rank.*/
    template <class Type, int rank> class rankType {};

    template <class Type> class rankType<Type, 0> {
    public:
        using type = Type;
    };

    template <class Type1, class Type2> class outerProduct {
    public:
        using type = typename rankType<typename feature<Type1>::elementType,
                                       int(feature<Type1>::rank + feature<Type2>::rank)>::type;
    };

    template <class Type1, class Type2> class crossProduct {
    public:
        using type = typename rankType<typename feature<Type1>::elementType,
                                       int(feature<Type1>::rank + feature<Type2>::rank - 1)>::type;
    };

    template <class Type1, class Type2> class innerProduct {
    public:
        using type = typename rankType<typename feature<Type1>::elementType,
                                       int(feature<Type1>::rank + feature<Type2>::rank - 2)>::type;
    };

} // namespace OpenHurricane