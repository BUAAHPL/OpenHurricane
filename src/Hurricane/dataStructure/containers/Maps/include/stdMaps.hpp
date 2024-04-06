/*!
 * \file stdMaps.hpp
 * \brief Header of stdMaps
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
#include "string.hpp"
#include <map>
#include <unordered_map>

namespace OpenHurricane {
    /**
     * \brief The template class of creating map.
     */
    template <typename keyType, typename valueType> class createMap {
    private:
        std::map<keyType, valueType> map_;

    public:
        /**\brief Constructors.*/
        createMap(const keyType &key, const valueType &val) { map_[key] = val; }

        createMap<keyType, valueType> &operator()(const keyType &key, const valueType &val) {
            map_[key] = val;
            return *this;
        }

        hur_nodiscard operator std::map<keyType, valueType>() { return map_; }
    };

    template <typename valueType>
    hur_nodiscard inline std::string stringMapDoc(const std::map<std::string, valueType> &strMap) {
        typename std::map<std::string, valueType>::const_iterator iter;
        std::string doc;
        doc = "\n";
        for (iter = strMap.begin(); iter != strMap.end(); ++iter) {
            doc += iter->first;
            doc += "\n";
        }
        return doc;
    }

    /**
     * \brief The template class of creating Unorderedmap.
     */
    template <typename keyType, typename valueType> class createUnorderedMap {
    private:
        std::unordered_map<keyType, valueType> unorderedMap_;

    public:
        /**\brief Constructors.*/
        createUnorderedMap(const keyType &key, const valueType &val) { unorderedMap_[key] = val; }

        createUnorderedMap<keyType, valueType> &operator()(const keyType &key,
                                                           const valueType &val) {
            unorderedMap_[key] = val;
            return *this;
        }

        hur_nodiscard operator std::unordered_map<keyType, valueType>() { return unorderedMap_; }
    };

    template <typename valueType>
    hur_nodiscard inline std::string
    stringMapDoc(const std::unordered_map<std::string, valueType> &strMap) {
        typename std::unordered_map<std::string, valueType>::const_iterator iter;
        std::string doc;
        doc = "\n";
        for (iter = strMap.begin(); iter != strMap.end(); ++iter) {
            doc += iter->first;
            doc += "\n";
        }
        return doc;
    }
} // namespace OpenHurricane