/*!
 * \file transfer.hpp
 * \brief Headers of the transfer.
 *        The subroutines and functions are in the <i>transfer.cpp</i> file.
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
#include "integerTransfer.hpp"
#include "realTransfer.hpp"
#include "tensorTransfer.hpp"
#include "transferTemplate.hpp"
#include "vectorTransfer.hpp"

namespace OpenHurricane {
    namespace fv {
        /*!\brief Transfer var between partition boundaries and between periodic boundaries.*/
        extern void transfer(cellRealArray &var, const bool onlyFirstLayer = false,
                             const bool block = false);

        /*!\brief Transfer var between partition boundaries.*/
        extern void cutTransfer(cellRealArray &var, const bool onlyFirstLayer = false,
                                const bool block = false);
        extern void cutTransfer(const runtimeMesh &mesh, realArray &var,
                                const bool onlyFirstLayer = false);
        extern void cutNonBlockingTransfer(const runtimeMesh &mesh, realArray &var,
                                           const bool onlyFirstLayer = false);

        /*!\brief Transfer var between periodic boundaries.*/
        extern void perTransfer(cellRealArray &var, const bool onlyFirstLayer = false,
                                const bool block = false);
        extern void perTransfer(const runtimeMesh &mesh, realArray &var,
                                const bool onlyFirstLayer = false);
        extern void perNonBlockingTransfer(const runtimeMesh &mesh, realArray &var,
                                           const bool onlyFirstLayer = false);

        /*!\brief Transfer var between partition boundaries and between periodic boundaries.*/
        extern void transfer(cellVectorArray &var, const bool onlyFirstLayer = false,
                             const bool block = false);

        /*!\brief Transfer var between partition boundaries.*/
        extern void cutTransfer(cellVectorArray &var, const bool onlyFirstLayer = false,
                                const bool block = false);
        extern void cutTransfer(const runtimeMesh &mesh, vectorArray &var,
                                const bool onlyFirstLayer = false);
        extern void cutNonBlockingTransfer(const runtimeMesh &mesh, vectorArray &var,
                                           const bool onlyFirstLayer = false);

        /*!\brief Transfer var between periodic boundaries.*/
        extern void perTransfer(cellVectorArray &var, const bool onlyFirstLayer = false,
                                const bool block = false);
        extern void perTransfer(const runtimeMesh &mesh, vectorArray &var,
                                const bool onlyFirstLayer = false);
        extern void perNonBlockingTransfer(const runtimeMesh &mesh, vectorArray &var,
                                           const bool onlyFirstLayer = false);

        /*!\brief Transfer var between partition boundaries and between periodic boundaries.*/
        extern void transfer(cellTensorArray &var, const bool onlyFirstLayer = false,
                             const bool block = false);
        
        /*!\brief Transfer var between partition boundaries.*/
        extern void cutTransfer(cellTensorArray &var, const bool onlyFirstLayer = false,
                                const bool block = false);
        extern void cutTransfer(const runtimeMesh &mesh, tensorArray &var,
                                const bool onlyFirstLayer = false);
        extern void cutNonBlockingTransfer(const runtimeMesh &mesh, tensorArray &var,
                                           const bool onlyFirstLayer = false);
        /*!\brief Transfer var between periodic boundaries.*/
        extern void perTransfer(cellTensorArray &var, const bool onlyFirstLayer = false,
                                const bool block = false);
        extern void perTransfer(const runtimeMesh &mesh, tensorArray &var,
                                const bool onlyFirstLayer = false);
        extern void perNonBlockingTransfer(const runtimeMesh &mesh, tensorArray &var,
                                           const bool onlyFirstLayer = false);
        //extern void transfer(cellSymmTensorArray& var);

        /*!\brief Transfer var between partition boundaries and between periodic boundaries.*/
        extern void transfer(const runtimeMesh &mesh, realArrayArray &var,
                             const bool onlyFirstLayer = false);
        extern void cutTransfer(const runtimeMesh &mesh, realArrayArray &var,
                                const bool onlyFirstLayer = false);
        extern void perTransfer(const runtimeMesh &mesh, realArrayArray &var,
                                const bool onlyFirstLayer = false);

    } // namespace fv
} // namespace OpenHurricane

#include "transfer.inl"