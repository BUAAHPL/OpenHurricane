/*!
 * \file uiListHashMap.hpp
 * \brief Header of unsigned integer List Hash Map.
 *       The subroutines and functions are in the <i>uiListHashMap.inl</i> file.
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
#include "HashMap.hpp"
#include "Lists.hpp"
#include "stdMaps.hpp"

namespace OpenHurricane {
    namespace HashFunctions {
        /**
         * \brief BKDR Hash Function.
         * This hash function comes from Brian Kernighan
         * and Dennis Ritchie's book "The C Programming Language".
         * A strange set of possible seeds which all constitute a pattern of 31,131,1313,13131,131313, etc.
         */
        template <class Type, uinteger seed = 131> class BKDRListHash {
        public:
            /**
             * \brief BKDR Hash Function.
             * This hash function comes from Brian Kernighan
             * and Dennis Ritchie's book "The C Programming Language".
             * A strange set of possible seeds which all constitute a pattern of 31,131,1313,13131,131313, etc.
             */
            hur_nodiscard inline uinteger operator()(const List<Type> &k) const {
                uinteger sum = 0;
                for (integer i = 0; i < k.size(); ++i) {
                    sum = sum * seed + k[i];
                }
                return sum;
            }
        };

        template <typename Type> class ListHashCompare {
        public:
            hur_nodiscard inline bool operator()(const List<Type> &lhs,
                                                 const List<Type> &rhs) const {
                return lhs == rhs;
            }
        };
    } // namespace HashFunctions

    /**
     * \brief List Hash map.
     */
    template <class keyElementType, class Type, uinteger seed = 131, bool usePrimeNumber = true>
    using ListHashMap =
        HashMap<List<keyElementType>, Type, HashFunctions::BKDRListHash<keyElementType, seed>,
                HashFunctions::ListHashCompare<keyElementType>, usePrimeNumber>;
    /**
     * \brief Integer List Hash map.
     */
    template <class Type, uinteger seed = 131, bool usePrimeNumber = true>
    using iListHashMap = ListHashMap<integer, Type, seed, usePrimeNumber>;

    /**
     * \brief Unsigned Integer List Hash map.
     */
    template <class Type, uinteger seed = 131, bool usePrimeNumber = true>
    using uiListHashMap = ListHashMap<uinteger, Type, seed, usePrimeNumber>;

    template <typename type, unsigned int seed, unsigned int modulus> class ListHash {
    public:
        /**
         * \brief BKDR Hash Function.
         * This hash function comes from Brian Kernighan
         * and Dennis Ritchie's book "The C Programming Language".
         * A strange set of possible seeds which all constitute a pattern of 31,131,1313,13131,131313, etc.
         */
        hur_nodiscard inline auto operator()(const List<type> &k) const {
            type sum = 0;
            for (integer i = 0; i < k.size(); ++i) {
                sum = sum * seed + k[i];
            }
            return std::hash<type>{}(sum % modulus);
        }
    };

    template <typename type> class ListHashCompare {
    public:
        hur_nodiscard inline bool operator()(const List<type> &lhs, const List<type> &rhs) const {
            return lhs == rhs;
        }
    };

    /**
     * \brief Integer List Hash table.
     */
    template <class Type, unsigned int seed = 131, unsigned int modulus = 100001651>
    using intListHashTable = std::unordered_map<integerList, Type, ListHash<integer, seed, modulus>,
                                                ListHashCompare<integer>>;

    /**
     * \brief Unsigned int List Hash table.
     */
    template <class Type, unsigned int seed = 131, unsigned int modulus = 100001651>
    using uintListHashTable =
        std::unordered_map<List<unsigned int>, Type, ListHash<unsigned int, seed, modulus>,
                           ListHashCompare<unsigned int>>;

    template <typename type, unsigned int modulus> class SetHash {
    public:
        hur_nodiscard inline auto operator()(const List<type> &k) const {
            type sum = 0;
            for (integer i = 0; i < k.size(); ++i) {
                sum += k[i];
            }
            if (modulus == 1) {
                return std::hash<type>{}(sum);
            } else {
                return std::hash<type>{}(sum % modulus);
            }
        }
    };

    template <typename type> class SetHashCompare {
    public:
        hur_nodiscard inline bool operator()(const List<type> &lhs, const List<type> &rhs) const {
            if (lhs.size() == rhs.size()) {
                bool equal = true;
                for (integer i = 0; i < lhs.size(); ++i) {
                    bool equal1 = false;
                    for (integer j = 0; j < rhs.size(); ++j) {
                        if (lhs[i] == rhs[j]) {
                            equal1 = true;
                            break;
                        }
                    }
                    if (!equal1) {
                        equal = false;
                        break;
                    }
                }
                return equal;
            }
            return false;
        }
    };

    /**
     * \brief Integer List Hash table.
     */
    template <class Type, unsigned int modulus = 1>
    using intSetHashTable =
        std::unordered_map<integerList, Type, SetHash<integer, modulus>, SetHashCompare<integer>>;

    /**
     * \brief Unsigned int List Hash table.
     */
    template <class Type, unsigned int modulus = 1>
    using uintSetHashTable =
        std::unordered_map<List<unsigned int>, Type, SetHash<unsigned int, modulus>,
                           SetHashCompare<unsigned int>>;

} // namespace OpenHurricane
