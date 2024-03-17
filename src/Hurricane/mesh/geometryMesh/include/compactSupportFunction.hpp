/*!
 * \file compactSupportFunction.hpp
 * \brief Header of compact support function for RBF.
 *       The subroutines and functions are in the <i>compactSupportFunction.cpp</i> file.
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

#include "OpenHurricane.hpp"

namespace OpenHurricane {
    /**
     * \brief The base class of compact support function for RBF.
     */
    class compactSupportFunction {
    protected:
        /*!\brief Support radius.*/
        real R_;

    public:
        declareClassNames;
        declareObjFty(compactSupportFunction, controller, (const controller &cont), (cont));

        /*!\brief Disallow null constructor.*/
        compactSupportFunction() = delete;

        /*!\brief Disallow copy constructor.*/
        compactSupportFunction(const compactSupportFunction &) = delete;

        /*!\brief Construct from controller.*/
        compactSupportFunction(const controller &cont);

        /*!\brief Disallow default bitwise assignment.*/
        void operator=(const compactSupportFunction &) = delete;

        static uniquePtr<compactSupportFunction> creator(const controller &cont);

        virtual ~compactSupportFunction() noexcept {}

    protected:
        hur_nodiscard virtual real f(const real psi) const = 0;

    public:
        hur_nodiscard inline real phi(const real r) const { return f(r / max(R_, tiny)); }
    };

    /**
     * \brief The base class of compact support function: CPC0 for RBF.
     */
    class CPC0 : public compactSupportFunction {
    public:
        declareClassNames;

        /*!\brief Disallow null constructor.*/
        CPC0() = delete;

        /*!\brief Disallow copy constructor.*/
        CPC0(const CPC0 &) = delete;

        /*!\brief Construct from controller.*/
        CPC0(const controller &cont);

        virtual ~CPC0() noexcept {}

    protected:
        hur_nodiscard virtual real f(const real psi) const;
    };

    /**
     * \brief The base class of compact support function: CPC2 for RBF.
     */
    class CPC2 : public compactSupportFunction {
    public:
        declareClassNames;

        /*!\brief Disallow null constructor.*/
        CPC2() = delete;

        /*!\brief Disallow copy constructor.*/
        CPC2(const CPC2 &) = delete;

        /*!\brief Construct from controller.*/
        CPC2(const controller &cont);

        virtual ~CPC2() noexcept {}

    protected:
        hur_nodiscard virtual real f(const real psi) const;
    };

    /**
     * \brief The base class of compact support function: CPC4 for RBF.
     */
    class CPC4 : public compactSupportFunction {
    public:
        declareClassNames;

        // Constructors

        /*!\brief Disallow null constructor.*/
        CPC4() = delete;

        /*!\brief Disallow copy constructor.*/
        CPC4(const CPC4 &) = delete;

        /*!\brief Construct from controller.*/
        CPC4(const controller &cont);

        virtual ~CPC4() noexcept {}

    protected:
        hur_nodiscard virtual real f(const real psi) const;
    };

    /**
     * \brief The base class of compact support function: CPC6 for RBF.
     */
    class CPC6 : public compactSupportFunction {
    public:
        declareClassNames;

        // Constructors

        /*!\brief Disallow null constructor.*/
        CPC6() = delete;

        /*!\brief Disallow copy constructor.*/
        CPC6(const CPC6 &) = delete;

        /*!\brief Construct from controller.*/
        CPC6(const controller &cont);

        virtual ~CPC6() noexcept {}

    protected:
        hur_nodiscard virtual real f(const real psi) const;
    };

    /**
     * \brief The base class of compact support function: CTPSC0 for RBF.
     */
    class CTPSC0 : public compactSupportFunction {
    public:
        declareClassNames;

        // Constructors

        /*!\brief Disallow null constructor.*/
        CTPSC0() = delete;

        /*!\brief Disallow copy constructor.*/
        CTPSC0(const CTPSC0 &) = delete;

        /*!\brief Construct from controller.*/
        CTPSC0(const controller &cont);

        virtual ~CTPSC0() noexcept {}

    protected:
        hur_nodiscard virtual real f(const real psi) const;
    };

    /**
     * \brief The base class of compact support function: CTPSC1 for RBF.
     */
    class CTPSC1 : public compactSupportFunction {
    public:
        declareClassNames;

        // Constructors

        /*!\brief Disallow null constructor.*/
        CTPSC1() = delete;

        /*!\brief Disallow copy constructor.*/
        CTPSC1(const CTPSC1 &) = delete;

        /*!\brief Construct from controller.*/
        CTPSC1(const controller &cont);

        virtual ~CTPSC1() noexcept {}

    protected:
        hur_nodiscard virtual real f(const real psi) const;
    };

    /**
     * \brief The base class of compact support function: CTPSC2a for RBF.
     */
    class CTPSC2a : public compactSupportFunction {
    public:
        declareClassNames;

        // Constructors

        /*!\brief Disallow null constructor.*/
        CTPSC2a() = delete;

        /*!\brief Disallow copy constructor.*/
        CTPSC2a(const CTPSC2a &) = delete;

        /*!\brief Construct from controller.*/
        CTPSC2a(const controller &cont);

        virtual ~CTPSC2a() noexcept {}

    protected:
        hur_nodiscard virtual real f(const real psi) const;
    };

    /**
     * \brief The base class of compact support function: CTPSC2b for RBF.
     */
    class CTPSC2b : public compactSupportFunction {
    public:
        declareClassNames;

        // Constructors

        /*!\brief Disallow null constructor.*/
        CTPSC2b() = delete;

        /*!\brief Disallow copy constructor.*/
        CTPSC2b(const CTPSC2b &) = delete;

        /*!\brief Construct from controller.*/
        CTPSC2b(const controller &cont);

        virtual ~CTPSC2b() noexcept {}

    protected:
        hur_nodiscard virtual real f(const real psi) const;
    };

} // namespace OpenHurricane
