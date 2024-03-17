/*!
 * \file objectFactory.hpp
 * \brief Header of object factory.
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

#include "errorAbort.hpp"
#include "setClassName.hpp"
#include "smartPointer.hpp"
#include <map>

/**
 * \brief Declare the object factory for the given base class: baseClass.
 * \param[in] baseClass - The base class.
 * \param[in] control - The scope of control.
 * \param[in] arguments - The arguments for the constructors of the given base class and its derived class.
 * \param[in] parmameters - The parmameters for the constructors of the given base class and its derived class.
 */
#define declareObjFty(baseClass, control, arguments, parmameters)                               \
    typedef uniquePtr<baseClass>(*control##creactorFuncPtr) arguments;                          \
    using control##creactorMapType = std::map<std::string, control##creactorFuncPtr>;           \
    static control##creactorMapType *control##creatorMapPtr_;                                   \
    static void control##createCreactorMap();                                                   \
    static void control##destroyCreactorMap();                                                  \
    class control##managementCreactorMap {                                                      \
    public:                                                                                     \
        control##managementCreactorMap() {                                                      \
            baseClass::control##createCreactorMap();                                            \
        }                                                                                       \
        virtual ~control##managementCreactorMap() noexcept {                                    \
            baseClass::control##destroyCreactorMap();                                           \
        }                                                                                       \
    };                                                                                          \
    template <class derivedClass>                                                               \
    class control##registerInCreactorMap : public control##managementCreactorMap {              \
    public:                                                                                     \
        static hur_nodiscard uniquePtr<baseClass> creatorFunction arguments {                   \
            return uniquePtr<baseClass>(new derivedClass parmameters);                          \
        }                                                                                       \
        control##registerInCreactorMap(const std::string &classType = derivedClass::className_) \
            : control##managementCreactorMap() {                                                \
            baseClass::control##creatorMapPtr_->emplace(classType, creatorFunction);            \
        }                                                                                       \
        virtual ~control##registerInCreactorMap() noexcept {}                                   \
    }

#define initializeCreatorMapPointer(baseClass, control) \
    typename baseClass::control##creactorMapType *baseClass::control##creatorMapPtr_ = nullptr;

#define defineCreateCreatorMap(baseClass, control)                  \
    void baseClass::control##createCreactorMap() {                  \
        if (control##creatorMapPtr_ == nullptr) {                   \
            control##creatorMapPtr_ = new control##creactorMapType; \
        }                                                           \
    }

#define defineDestroyCreactorMap(baseClass, control) \
    void baseClass::control##destroyCreactorMap() {  \
        HurDelete(control##creatorMapPtr_);          \
    }

/**
 * \brief Create the object factory for the given base class: baseClass.
 * \param[in] baseClass - The base class.
 * \param[in] control - The scope of control.
 */
#define createObjFty(baseClass, control)             \
    initializeCreatorMapPointer(baseClass, control); \
    defineCreateCreatorMap(baseClass, control);      \
    defineDestroyCreactorMap(baseClass, control);

/**
 * \brief Create the object factory for the given base class: baseClass.
 * \param[in] baseClass - The base class.
 * \param[in] control - The scope of control.
 */
#define createObjFtyTmpl(baseClass, control)                     \
    template <> initializeCreatorMapPointer(baseClass, control); \
    template <> defineCreateCreatorMap(baseClass, control);      \
    template <> defineDestroyCreactorMap(baseClass, control);

/**
 * \brief Register the derived class into the object factory of the given base class.
 * \param[in] baseClass - The base class.
 * \param[in] derivedClass - The derived class.
 * \param[in] control - The scope of control.
 */
#define registerObjFty(baseClass, derivedClass, control)    \
    baseClass::control##registerInCreactorMap<derivedClass> \
        baseClass##control##derivedClass##RegisterInMap

#define defineInObjCreator(baseClass, derivedClassName, control, parmameters) \
    control##creactorMapType::iterator control##iter =                        \
        control##creatorMapPtr_->find(derivedClassName);                      \
    if (control##iter == control##creatorMapPtr_->end()) {                    \
        std::string errMsg;                                                   \
        errMsg = "Unknown ";                                                  \
        errMsg += #baseClass;                                                 \
        errMsg += "type: ";                                                   \
        errMsg += derivedClassName;                                           \
        errMsg += "\nValid ";                                                 \
        errMsg += #baseClass;                                                 \
        errMsg += "types of current program are: ";                           \
        errMsg += stringMapDoc(*control##creatorMapPtr_);                     \
        errorAbortStr(errMsg);                                                \
    }                                                                         \
    return (control##iter->second)parmameters

#define defineInObjCreatorTmpl(baseClass, derivedClassName, control, parmameters) \
    typename control##creactorMapType::iterator control##iter =                   \
        control##creatorMapPtr_->find(derivedClassName);                          \
    if (control##iter == control##creatorMapPtr_->end()) {                        \
        std::string errMsg;                                                       \
        errMsg = "Unknown ";                                                      \
        errMsg += #baseClass;                                                     \
        errMsg += "type: ";                                                       \
        errMsg += derivedClassName;                                               \
        errMsg += "\nValid ";                                                     \
        errMsg += #baseClass;                                                     \
        errMsg += "types of current program are: ";                               \
        errMsg += stringMapDoc(*control##creatorMapPtr_);                         \
        errorAbortStr(errMsg);                                                    \
    }                                                                             \
    return (control##iter->second)parmameters
