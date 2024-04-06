/*!
 * \file transferTemplate.hpp
 * \brief Headers of the transfer template.
 *        The subroutines and functions are in the <i>transferTemplate.inl</i> file.
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
#include "fVArraysInclude.hpp"

namespace OpenHurricane {
    template <class Type> class processTransfer {
    private:
        const runtimeMesh &mesh_;

        Array<Type> &var_;

        bool isBlock_;

        bool onlyFirstLayer_;

        List<HurMPI::Request> requestCut_;
        PtrList<cutProcessReceiv<Type>> cutPRecvs_;
        PtrList<cutProcessSend<Type>> cutPSends_;

        List<HurMPI::Request> requestPer_;
        PtrList<perProcessReceiv<Type>> perPRecvs_;
        PtrList<perProcessSend<Type>> perPSends_;

    public:
        // Constructors

        /**
         * \brief Disallow null constructor.
         */
        processTransfer() = delete;

        /**
         * \brief Construct from.
         */
        inline processTransfer(const runtimeMesh &mesh, Array<Type> &var, const bool isBlock = true,
                               const bool onlyFirstLayer = false);

        /**
         * \brief Disallow copy constructor.
         */
        processTransfer(const processTransfer &) = delete;

        /**
         * \brief Destructor.
         */
        inline ~processTransfer() noexcept;

        void clear() noexcept;

        void transferInit();

        void transferInit(const integer layerI);

        void transferring();

        void transferring(const integer layerI);

        inline void waitAll();

    private:
        void transferringBlock();
        void transferringBlock(const integer layerI);
        void cutTransferringBlock(const cutZone &cuts);
        void perTransferringBlock(const perZone &pers);

    private:
        void cutNonBlockingTransferInit();
        void cutNonBlockingTransferInit(const integer layerI);

        void perNonBlockingTransferInit();
        void perNonBlockingTransferInit(const integer layerI);

        inline void cutNonBlockingTransferring();
        inline void perNonBlockingTransferring();

    public:
        /**
         * \brief Disallow bitwise assignment.
         */
        void operator=(const processTransfer &) = delete;
    };
} // namespace OpenHurricane

#include "transferTemplate.inl"