/*!
 * \file creatingFunctionMap.hpp
 * \brief Header of creating function map.
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

#include "stdMaps.hpp"
#include "string.hpp"

#define creatingFunctionMap(newClass, keyType, returnType, argList, errFunc)          \
    typedef returnType(*newClass##FuncMapPtr) argList;                                \
                                                                                      \
    class newClass {                                                                  \
    private:                                                                          \
        std::map<keyType, newClass##FuncMapPtr> funcMap_;                             \
                                                                                      \
    public:                                                                           \
        inline newClass() : funcMap_() {}                                             \
                                                                                      \
        newClass(const newClass &) = delete;                                          \
                                                                                      \
        inline ~newClass() noexcept {}                                                \
                                                                                      \
        inline bool found(const keyType &key) const {                                 \
            if (funcMap_.empty()) {                                                   \
                return false;                                                         \
            }                                                                         \
            const auto iter = funcMap_.find(key);                                     \
            if (iter != funcMap_.end()) {                                             \
                return true;                                                          \
            }                                                                         \
            return false;                                                             \
        }                                                                             \
                                                                                      \
        inline void clear() noexcept {                                                \
            funcMap_.clear();                                                         \
        }                                                                             \
                                                                                      \
        inline newClass##FuncMapPtr find(const keyType &key) const {                  \
            const auto iter = funcMap_.find(key);                                     \
            return iter->second;                                                      \
        }                                                                             \
                                                                                      \
        inline void addFunc(const keyType &key, newClass##FuncMapPtr newFunc) {       \
            if (found(key)) {                                                         \
                errFunc(("The function of " + key + " is already added to the map")); \
            }                                                                         \
            funcMap_.emplace(key, newFunc);                                           \
        }                                                                             \
    };
