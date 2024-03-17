/*!
 * \file Lists.hpp
 * \brief Header of lists
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
#include "List.hpp"
#include "real.hpp"
#include "bools.hpp"
#include "string.hpp"
#include "integer.hpp"
#include "sinteger.hpp"
#include "vectors.hpp"
#include "tensors.hpp"

namespace OpenHurricane {
    
    using realList = List<real>;
    using realListList = List<realList>;

    using boolList = List<bool>;
    using boolListList = List<List<bool>>;    

    using integerList = List<integer>;
    using integerListList = List<integerList>;
    using integerListListList = List<integerListList>;

    using sintegerList = List<sinteger>;
    using sintegerListList = List<sintegerList>;
    using sintegerListListList = List<sintegerListList>;
        
    using vectorList = List<vector>;
    using integerVectorList = List<integerVector>;

    using vector2DList = List<vector2D>;
    using integerVector2DList = List<integerVector2D>;
    
    using tensorList = List<tensor>;
    using symmTensorList = List<symmTensor>;
    using sphericalTensorList = List<sphericalTensor>;

    using stdStringList = List<std::string>;
    using stringList = List<string>;

    hur_nodiscard stdStringList stringToLineList(const std::string &str);

    hur_nodiscard stringList stringToLineList(const string &str);

    /**\brief Split string to string list by token given by character c.*/
    void split(const std::string &str, stdStringList &taken, const char *c = " ");

    void split(const std::string &str, stdStringList &taken, const std::string &c = " ");

    void split(const std::string &str, stringList &taken, const char *c = " ");

    void split(const std::string &str, stringList &taken, const std::string &c = " ");

} // namespace OpenHurricane