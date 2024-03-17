/*!
 * \file preset.hpp
 * \brief Headers of pre-set code.
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

#ifdef _MSC_VER // for MSVC
#define hur_forceinline __forceinline
#define hur_inline inline
#define hur_restrict __restrict
#elif defined(__GNUC__) || (defined(_ICC) && (__ICC >= 600)) // for gcc on Linux/Apple OS X
#define hur_forceinline __inline__ __attribute__((always_inline))
#define hur_inline inline
#define hur_restrict __restrict__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600))
#define hur_forceinline inline
#define hur_inline inline
#define hur_restrict __restrict
#else
#define hur_forceinline inline
#define hur_inline inline
#define hur_restrict __restrict
#endif

#if defined(_MSC_VER) && (_MSC_VER >= 1800)
#define HUR_FUNCTION __FUNCTION__
#define HUR_FUNC __func__
#elif defined(__GNUC__) || (defined(_ICC) && (__ICC >= 600))
#define HUR_FUNCTION __PRETTY_FUNCTION__
#define HUR_FUNC __FUNCTION__
#elif defined(__FUNCSIG__)
#define HUR_FUNCTION __FUNCSIG__
#define HUR_FUNC HUR_FUNCTION
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600))
#define HUR_FUNCTION __FUNCTION__
#define HUR_FUNC HUR_FUNCTION
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
#define HUR_FUNCTION __func__
#define HUR_FUNC HUR_FUNCTION
#elif defined(__cplusplus) && (__cplusplus >= 201103)
#define HUR_FUNCTION __func__
#else
#define HUR_FUNCTION "(unknown)"
#endif

#if defined(__GNUC__)
#define hur_noreturn __attribute__((noreturn))
#define hur_deprecated __attribute__((deprecated))
#define hur_nodiscard
#if __GNUC__ >= 6
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#if __GNUC__ == 7
// See: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=89325
#pragma GCC diagnostic ignored "-Wattributes"
#endif

#elif defined(_MSC_VER)
#define hur_noreturn __declspec(noreturn)
#define hur_deprecated __declspec(deprecated)
#if (_MSC_VER >= 1930)
#define hur_nodiscard [[nodiscard]]
#else
#define hur_nodiscard [[nodiscard]]
#endif
#pragma warning(disable : 26812)

#else
#define hur_noreturn
#define hur_deprecated
#define hur_nodiscard
#endif

#define HUR_LINE __LINE__

#define HUR_FILE __FILE__

#ifndef HurDelete
/** \brief Delete the pointer if the pointer is not null. */
#define HurDelete(ptr)    \
    if (ptr != nullptr) { \
        delete ptr;       \
        ptr = nullptr;    \
    }

#endif // !HurDelete

#ifndef HurDeleteDynArray
/** \brief Delete the pointer if the pointer is not null. */
#define HurDeleteDynArray(ptr) \
    if (ptr != nullptr) {      \
        delete[] ptr;          \
        ptr = nullptr;         \
    }

#endif // !HurDeleteDynArray

#ifndef HUR_HAS_CXX17
#   ifdef _MSC_VER
#       if _MSVC_LANG > 201402L
#           define HUR_HAS_CXX17 1
#       else
#           define HUR_HAS_CXX17 0
#       endif
#   else
#       if __cplusplus > 201402L
#           define HUR_HAS_CXX17 1
#       else 
#           define HUR_HAS_CXX17 0
#       endif
#   endif
#endif // _HAS_CXX17

#if HUR_HAS_CXX17
#define hur_inline_var inline
#else
#define hur_inline_var
#endif
