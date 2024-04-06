#include "transfer.hpp"
/*!
 * \file transfer.inl
 * \brief In-Line subroutines of the <i>transfer.hpp</i> file.
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

inline void OpenHurricane::fv::transfer(cellRealArray &var,
                                        const bool onlyFirstLayer,
                                        const bool block) {
    if (HurMPI::parRun()) {
        //20210408 ��˼�� �������߽�����ڱ߽紫�����ݷֿ�
        cutTransfer(var, onlyFirstLayer, block);
    }
    perTransfer(var, onlyFirstLayer, block);
}

inline void OpenHurricane::fv::cutTransfer(cellRealArray &var,
                                           const bool onlyFirstLayer,
                                           const bool block) {
    if (block) {
        cutTransfer(var.mesh(), var, onlyFirstLayer);
    } else {
        cutNonBlockingTransfer(var.mesh(), var, onlyFirstLayer);
    }
}

inline void OpenHurricane::fv::perTransfer(cellRealArray &var,
                                           const bool onlyFirstLayer,
                                           const bool block) {
    if (block) {
        perTransfer(var.mesh(), var, onlyFirstLayer);
    } else {
        perNonBlockingTransfer(var.mesh(), var, onlyFirstLayer);
    }
}

inline void OpenHurricane::fv::transfer(cellVectorArray &var,
                                        const bool onlyFirstLayer,
                                        const bool block) {
    //20210408 ��˼�� �������߽�����ڱ߽紫�����ݷֿ�
    if (HurMPI::parRun()) {
        cutTransfer(var, onlyFirstLayer, block);
    }
    perTransfer(var, onlyFirstLayer, block);
}

inline void OpenHurricane::fv::cutTransfer(cellVectorArray &var,
                                           const bool onlyFirstLayer,
                                           const bool block) {
    if (block) {
        cutTransfer(var.mesh(), var, onlyFirstLayer);
    } else {
        cutNonBlockingTransfer(var.mesh(), var, onlyFirstLayer);
    }
}

inline void OpenHurricane::fv::perTransfer(cellVectorArray &var,
                                           const bool onlyFirstLayer,
                                           const bool block) {
    if (block) {
        perTransfer(var.mesh(), var, onlyFirstLayer);
    } else {
        perNonBlockingTransfer(var.mesh(), var, onlyFirstLayer);
    }
}

inline void OpenHurricane::fv::transfer(cellTensorArray &var,
                                        const bool onlyFirstLayer,
                                        const bool block) {
    //20210408 ��˼�� �������߽�����ڱ߽紫�����ݷֿ�
    if (HurMPI::parRun()) {
        cutTransfer(var, onlyFirstLayer, block);
    }
    perTransfer(var, onlyFirstLayer, block);
}

inline void OpenHurricane::fv::cutTransfer(cellTensorArray &var,
                                           const bool onlyFirstLayer,
                                           const bool block) {
    if (block) {
        cutTransfer(var.mesh(), var, onlyFirstLayer);
    } else {
        cutNonBlockingTransfer(var.mesh(), var, onlyFirstLayer);
    }
}

inline void OpenHurricane::fv::perTransfer(cellTensorArray &var,
                                           const bool onlyFirstLayer,
                                           const bool block) {
    if (block) {
        perTransfer(var.mesh(), var, onlyFirstLayer);
    } else {
        perNonBlockingTransfer(var.mesh(), var, onlyFirstLayer);
    }
}

inline void OpenHurricane::fv::transfer(const runtimeMesh &mesh,
                                        realArrayArray &var,
                                        const bool onlyFirstLayer) {
    cutTransfer(mesh, var, onlyFirstLayer);
    perTransfer(mesh, var, onlyFirstLayer);
}