/*!
 * \file transfer.cpp
 * \brief Main subroutines for transfer.
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
#include "transfer.hpp"

//
//void OpenHurricane::fv::transfer(cellRealArray& var, const bool onlyFirstLayer)
//{
//	//20210408 ��˼�� �������߽�����ڱ߽紫�����ݷֿ�
//	cutTransfer(var, onlyFirstLayer);
//	perTransfer(var, onlyFirstLayer);
//	/*const auto& mesh = var.mesh();
//	const integer ghostLayer = mesh.ghostCellsType();
//	const auto& cutZones_ = mesh.cutZones();
//	const auto& perZones_ = mesh.perZones();
//	integer czs = cutZones_.size() / ghostLayer;
//	if (czs != 0)
//	{
//		for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++)
//		{
//			if (onlyFirstLayer && layerI >= 1)
//			{
//				break;
//			}
//			for (integer i = czs * layerI; i < (layerI + 1) * czs; i++)
//			{
//				cutZones_[i].transfer<Array, real>(var);
//			}
//		}
//	}
//	czs = perZones_.size() / ghostLayer;
//	if (czs != 0)
//	{
//		for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++)
//		{
//			if (onlyFirstLayer && layerI >= 1)
//			{
//				break;
//			}
//			for (integer i = czs * layerI; i < (layerI + 1) * czs; i++)
//			{
//				perZones_[i].transfer<Array, real>(var);
//			}
//		}
//	}*/
//}

//void OpenHurricane::fv::cutTransfer(cellRealArray& var, const bool onlyFirstLayer)
//{
//	const auto& mesh = var.mesh();
//	const integer ghostLayer = mesh.ghostCellsType();
//	const auto& cutZones_ = mesh.cutZones();
//	integer czs = cutZones_.size() / ghostLayer;
//	if (czs != 0)
//	{
//		for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++)
//		{
//			if (onlyFirstLayer && layerI >= 1)
//			{
//				break;
//			}
//			for (integer i = czs * layerI; i < (layerI + 1) * czs; i++)
//			{
//				cutZones_[i].transfer<Array, real>(var);
//			}
//		}
//	}
//}

void OpenHurricane::fv::cutTransfer(const runtimeMesh &mesh, realArray &var,
                                const bool onlyFirstLayer) {
    if (!HurMPI::parRun()) {
        return;
    }
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the total cells (including "
                   "ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &cutZones_ = mesh.cutZones();
    integer czs = cutZones_.size() / ghostLayer;
    if (czs != 0) {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            if (onlyFirstLayer && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                cutZones_[i].transfer<Array, real>(var);
            }
        }
    } else if (cutZones_.size() != 0) {
        errorAbortStr(("The size of cutZone is " + toString(cutZones_.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
}

void OpenHurricane::fv::cutNonBlockingTransfer(const runtimeMesh &mesh, realArray &var,
                                           const bool onlyFirstLayer) {
    if (!HurMPI::parRun()) {
        return;
    }
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const auto &cutZones = mesh.cutZones();
    const auto &processCut = mesh.processCut();
    const integer ghostLayer = mesh.ghostCellsType();
    integer czs = cutZones.size() / ghostLayer;
    if ((czs == 0 || ghostLayer == 0) && cutZones.size() != 0) {
        errorAbortStr(("The size of cutZone is " + toString(cutZones.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
    integer sizeSend = Zero;
    integer sizeRecv = Zero;
    if (onlyFirstLayer) {
        sizeSend += processCut.zoneIdSend()[0].size();
        sizeRecv += processCut.zoneIdRecv()[0].size();
    } else {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            sizeSend += processCut.zoneIdSend()[layerI].size();
            sizeRecv += processCut.zoneIdRecv()[layerI].size();
        }
    }

    List<HurMPI::Request> requests(sizeSend + sizeRecv);

    PtrList<cutProcessReceiv<real>> cutPRecvs(sizeRecv);
    PtrList<cutProcessSend<real>> cutPSends(sizeSend);

    integer cnt = 0;
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutRecv = processCut.zoneIdRecv()[layerI];
        for (integer i = 0; i < cutRecv.size(); ++i) {
            const auto zid = cutRecv[i];
            cutPRecvs.set(i, new cutProcessReceiv<real>(cutZones[zid]));
            cutPRecvs[i].recvNonBlock(&requests[cnt++]);
        }
    }
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutSend = processCut.zoneIdSend()[layerI];
        for (integer i = 0; i < cutSend.size(); ++i) {
            const auto zid = cutSend[i];
            cutPSends.set(i, new cutProcessSend<real>(cutZones[zid]));
            cutPSends[i].writeBUf(var);
            cutPSends[i].sendNonBlock(&requests[cnt++]);
        }
    }
    HurMPI::waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    for (integer i = 0; i < cutPRecvs.size(); ++i) {
        cutPRecvs[i].readBUf(var);
    }
}

void OpenHurricane::fv::perTransfer(const runtimeMesh &mesh, realArray &var,
                                const bool onlyFirstLayer) {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &perZones_ = mesh.perZones();
    integer czs = perZones_.size() / ghostLayer;
    if (czs != 0) {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            if (onlyFirstLayer && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                perZones_[i].transfer<Array, real>(var);
            }
        }
    } else if (perZones_.size() != 0) {
        errorAbortStr(("The size of periodic zones is " + toString(perZones_.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
}

void OpenHurricane::fv::perNonBlockingTransfer(const runtimeMesh &mesh, realArray &var,
                                           const bool onlyFirstLayer) {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &perZones = mesh.perZones();
    if (perZones.size() == 0) {
        return;
    }
    integer czs = perZones.size() / ghostLayer;
    if ((czs == 0 || ghostLayer == 0) && perZones.size() != 0) {
        errorAbortStr(("The size of cutZone is " + toString(perZones.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }

    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
            if (perZones[i].sendProc() == perZones[i].receivProc()) {
                perZones[i].transfer<Array, real>(var);
            }
        }
    }

    if (!HurMPI::parRun()) {
        return;
    }

    const auto &processPer = mesh.processPer();
    integer sizeSend = Zero;
    integer sizeRecv = Zero;
    if (onlyFirstLayer) {
        sizeSend += processPer.zoneIdSend()[0].size();
        sizeRecv += processPer.zoneIdRecv()[0].size();
    } else {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            sizeSend += processPer.zoneIdSend()[layerI].size();
            sizeRecv += processPer.zoneIdRecv()[layerI].size();
        }
    }

    List<HurMPI::Request> requests(sizeSend + sizeRecv);

    PtrList<perProcessReceiv<real>> perPRecvs(sizeRecv);
    PtrList<perProcessSend<real>> perPSends(sizeSend);

    integer cnt = 0;
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutRecv = processPer.zoneIdRecv()[layerI];
        for (integer i = 0; i < cutRecv.size(); ++i) {
            const auto zid = cutRecv[i];
            perPRecvs.set(i, new perProcessReceiv<real>(perZones[zid]));
            perPRecvs[i].recvNonBlock(&requests[cnt++]);
        }
    }
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutSend = processPer.zoneIdSend()[layerI];
        for (integer i = 0; i < cutSend.size(); ++i) {
            const auto zid = cutSend[i];
            perPSends.set(i, new perProcessSend<real>(perZones[zid]));
            perPSends[i].writeBUf(var);
            perPSends[i].sendNonBlock(&requests[cnt++]);
        }
    }
    HurMPI::waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    for (integer i = 0; i < perPRecvs.size(); ++i) {
        perPRecvs[i].readBUf(var);
    }
}

void OpenHurricane::fv::cutTransfer(const runtimeMesh &mesh, vectorArray &var,
                                const bool onlyFirstLayer) {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &cutZones_ = mesh.cutZones();
    integer czs = cutZones_.size() / ghostLayer;
    if (czs != 0) {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            if (onlyFirstLayer && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                cutZones_[i].transferVS<Array, vector>(var);
            }
        }
    } else if (cutZones_.size() != 0) {
        errorAbortStr(("The size of cut zones is " + toString(cutZones_.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
}

void OpenHurricane::fv::cutNonBlockingTransfer(const runtimeMesh &mesh, vectorArray &var,
                                           const bool onlyFirstLayer) {
    if (!HurMPI::parRun()) {
        return;
    }
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const auto &cutZones = mesh.cutZones();
    const auto &processCut = mesh.processCut();
    const integer ghostLayer = mesh.ghostCellsType();
    integer czs = cutZones.size() / ghostLayer;
    if ((czs == 0 || ghostLayer == 0) && cutZones.size() != 0) {
        errorAbortStr(("The size of cutZone is " + toString(cutZones.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
    integer sizeSend = Zero;
    integer sizeRecv = Zero;
    if (onlyFirstLayer) {
        sizeSend += processCut.zoneIdSend()[0].size();
        sizeRecv += processCut.zoneIdRecv()[0].size();
    } else {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            sizeSend += processCut.zoneIdSend()[layerI].size();
            sizeRecv += processCut.zoneIdRecv()[layerI].size();
        }
    }

    List<HurMPI::Request> requests(sizeSend + sizeRecv);

    PtrList<cutProcessReceiv<vector>> cutPRecvs(sizeRecv);
    PtrList<cutProcessSend<vector>> cutPSends(sizeSend);

    integer cnt = 0;
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutRecv = processCut.zoneIdRecv()[layerI];
        for (integer i = 0; i < cutRecv.size(); ++i) {
            const auto zid = cutRecv[i];
            cutPRecvs.set(i, new cutProcessReceiv<vector>(cutZones[zid]));
            cutPRecvs[i].recvNonBlock(&requests[cnt++]);
        }
    }
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutSend = processCut.zoneIdSend()[layerI];
        for (integer i = 0; i < cutSend.size(); ++i) {
            const auto zid = cutSend[i];
            cutPSends.set(i, new cutProcessSend<vector>(cutZones[zid]));
            cutPSends[i].writeBUf(var);
            cutPSends[i].sendNonBlock(&requests[cnt++]);
        }
    }
    HurMPI::waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    for (integer i = 0; i < cutPRecvs.size(); ++i) {
        cutPRecvs[i].readBUf(var);
    }
}

void OpenHurricane::fv::perTransfer(const runtimeMesh &mesh, vectorArray &var,
                                const bool onlyFirstLayer) {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &perZones_ = mesh.perZones();
    integer czs = perZones_.size() / ghostLayer;
    if (czs != 0) {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            if (onlyFirstLayer && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                perZones_[i].transferVS<Array, vector>(var);
            }
        }
    } else if (perZones_.size() != 0) {
        errorAbortStr(("The size of periodic zones is " + toString(perZones_.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
}

void OpenHurricane::fv::perNonBlockingTransfer(const runtimeMesh &mesh, vectorArray &var,
                                           const bool onlyFirstLayer) {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &perZones = mesh.perZones();
    integer czs = perZones.size() / ghostLayer;
    if (czs != 0) {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            if (onlyFirstLayer && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                if (perZones[i].sendProc() == perZones[i].receivProc()) {
                    perZones[i].transferVS(var);
                }
            }
        }
    } else if (perZones.size() != 0) {
        errorAbortStr(("The size of periodic zones is " + toString(perZones.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }

    if (!HurMPI::parRun()) {
        return;
    }

    const auto &processPer = mesh.processPer();
    integer sizeSend = Zero;
    integer sizeRecv = Zero;
    if (onlyFirstLayer) {
        sizeSend += processPer.zoneIdSend()[0].size();
        sizeRecv += processPer.zoneIdRecv()[0].size();
    } else {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            sizeSend += processPer.zoneIdSend()[layerI].size();
            sizeRecv += processPer.zoneIdRecv()[layerI].size();
        }
    }

    List<HurMPI::Request> requests(sizeSend + sizeRecv);

    PtrList<perProcessReceiv<vector>> perPRecvs(sizeRecv);
    PtrList<perProcessSend<vector>> perPSends(sizeSend);

    integer cnt = 0;
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutRecv = processPer.zoneIdRecv()[layerI];
        for (integer i = 0; i < cutRecv.size(); ++i) {
            const auto zid = cutRecv[i];
            perPRecvs.set(i, new perProcessReceiv<vector>(perZones[zid]));
            perPRecvs[i].recvNonBlock(&requests[cnt++]);
        }
    }
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutSend = processPer.zoneIdSend()[layerI];
        for (integer i = 0; i < cutSend.size(); ++i) {
            const auto zid = cutSend[i];
            perPSends.set(i, new perProcessSend<vector>(perZones[zid]));
            perPSends[i].writeBUf(var);
            perPSends[i].sendNonBlock(&requests[cnt++]);
        }
    }
    HurMPI::waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    for (integer i = 0; i < perPRecvs.size(); ++i) {
        perPRecvs[i].readBUf(var);
    }
}

void OpenHurricane::fv::cutTransfer(const runtimeMesh &mesh, tensorArray &var,
                                const bool onlyFirstLayer) {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &cutZones_ = mesh.cutZones();
    integer czs = cutZones_.size() / ghostLayer;
    if (czs != 0) {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            if (onlyFirstLayer && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                cutZones_[i].transferVS<Array, tensor>(var);
            }
        }
    } else if (cutZones_.size() != 0) {
        errorAbortStr(("The size of cut zones is " + toString(cutZones_.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
}

void OpenHurricane::fv::cutNonBlockingTransfer(const runtimeMesh &mesh, tensorArray &var,
                                           const bool onlyFirstLayer) {
    if (!HurMPI::parRun()) {
        return;
    }
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const auto &cutZones = mesh.cutZones();
    const auto &processCut = mesh.processCut();
    const integer ghostLayer = mesh.ghostCellsType();
    integer czs = cutZones.size() / ghostLayer;
    if ((czs == 0 || ghostLayer == 0) && cutZones.size() != 0) {
        errorAbortStr(("The size of cutZone is " + toString(cutZones.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
    integer sizeSend = Zero;
    integer sizeRecv = Zero;
    if (onlyFirstLayer) {
        sizeSend += processCut.zoneIdSend()[0].size();
        sizeRecv += processCut.zoneIdRecv()[0].size();
    } else {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            sizeSend += processCut.zoneIdSend()[layerI].size();
            sizeRecv += processCut.zoneIdRecv()[layerI].size();
        }
    }

    List<HurMPI::Request> requests(sizeSend + sizeRecv);

    PtrList<cutProcessReceiv<tensor>> cutPRecvs(sizeRecv);
    PtrList<cutProcessSend<tensor>> cutPSends(sizeSend);

    integer cnt = 0;
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutRecv = processCut.zoneIdRecv()[layerI];
        for (integer i = 0; i < cutRecv.size(); ++i) {
            const auto zid = cutRecv[i];
            cutPRecvs.set(i, new cutProcessReceiv<tensor>(cutZones[zid]));
            cutPRecvs[i].recvNonBlock(&requests[cnt++]);
        }
    }
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutSend = processCut.zoneIdSend()[layerI];
        for (integer i = 0; i < cutSend.size(); ++i) {
            const auto zid = cutSend[i];
            cutPSends.set(i, new cutProcessSend<tensor>(cutZones[zid]));
            cutPSends[i].writeBUf(var);
            cutPSends[i].sendNonBlock(&requests[cnt++]);
        }
    }
    HurMPI::waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    for (integer i = 0; i < cutPRecvs.size(); ++i) {
        cutPRecvs[i].readBUf(var);
    }
}

void OpenHurricane::fv::perTransfer(const runtimeMesh &mesh, tensorArray &var,
                                const bool onlyFirstLayer) {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &perZones_ = mesh.perZones();
    integer czs = perZones_.size() / ghostLayer;
    if (czs != 0) {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            if (onlyFirstLayer && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                perZones_[i].transferTensorTran(var);
            }
        }
    } else if (perZones_.size() != 0) {
        errorAbortStr(("The size of periodic zones is " + toString(perZones_.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
}

void OpenHurricane::fv::perNonBlockingTransfer(const runtimeMesh &mesh, tensorArray &var,
                                           const bool onlyFirstLayer) {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &perZones = mesh.perZones();
    integer czs = perZones.size() / ghostLayer;
    if (czs != 0) {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            if (onlyFirstLayer && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                if (perZones[i].sendProc() == perZones[i].receivProc()) {
                    perZones[i].transferTensorTran(var);
                }
            }
        }
    } else if (perZones.size() != 0) {
        errorAbortStr(("The size of periodic zones is " + toString(perZones.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }

    if (!HurMPI::parRun()) {
        return;
    }

    const auto &processPer = mesh.processPer();
    integer sizeSend = Zero;
    integer sizeRecv = Zero;
    if (onlyFirstLayer) {
        sizeSend += processPer.zoneIdSend()[0].size();
        sizeRecv += processPer.zoneIdRecv()[0].size();
    } else {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            sizeSend += processPer.zoneIdSend()[layerI].size();
            sizeRecv += processPer.zoneIdRecv()[layerI].size();
        }
    }

    List<HurMPI::Request> requests(sizeSend + sizeRecv);

    PtrList<perProcessReceiv<tensor>> perPRecvs(sizeRecv);
    PtrList<perProcessSend<tensor>> perPSends(sizeSend);

    integer cnt = 0;
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutRecv = processPer.zoneIdRecv()[layerI];
        for (integer i = 0; i < cutRecv.size(); ++i) {
            const auto zid = cutRecv[i];
            perPRecvs.set(i, new perProcessReceiv<tensor>(perZones[zid]));
            perPRecvs[i].recvNonBlock(&requests[cnt++]);
        }
    }
    for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
        if (onlyFirstLayer && layerI >= 1) {
            break;
        }
        const auto &cutSend = processPer.zoneIdSend()[layerI];
        for (integer i = 0; i < cutSend.size(); ++i) {
            const auto zid = cutSend[i];
            perPSends.set(i, new perProcessSend<tensor>(perZones[zid]));
            perPSends[i].writeBUf(var);
            perPSends[i].sendNonBlock(&requests[cnt++]);
        }
    }
    HurMPI::waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    for (integer i = 0; i < perPRecvs.size(); ++i) {
        perPRecvs[i].readBUf(var);
    }
}

void OpenHurricane::fv::cutTransfer(const runtimeMesh &mesh, realArrayArray &var,
                                const bool onlyFirstLayer) {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &cutZones_ = mesh.cutZones();
    integer czs = cutZones_.size() / ghostLayer;
    if (czs != 0) {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            if (onlyFirstLayer && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                cutZones_[i].transfer<Array, realArray>(var);
            }
        }
    } else if (cutZones_.size() != 0) {
        errorAbortStr(("The size of cut zones is " + toString(cutZones_.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
}

void OpenHurricane::fv::perTransfer(const runtimeMesh &mesh, realArrayArray &var,
                                const bool onlyFirstLayer) {
    if (var.size() < mesh.nTotalCells()) {
        LFatal("Cannot transfer array of which size is less than the "
                   "total cells (including ghost cells)");
        return;
    }
    const integer ghostLayer = mesh.ghostCellsType();
    const auto &perZones_ = mesh.perZones();
    integer czs = perZones_.size() / ghostLayer;
    if (czs != 0) {
        for (integer layerI = 0; layerI < mesh.ghostCellsType(); layerI++) {
            if (onlyFirstLayer && layerI >= 1) {
                break;
            }
            for (integer i = czs * layerI; i < (layerI + 1) * czs; i++) {
                perZones_[i].transfer<Array, realArray>(var);
            }
        }
    } else if (perZones_.size() != 0) {
        errorAbortStr(("The size of periodic zones is " + toString(perZones_.size()) +
                       ", but the layer of ghost cells is " + toString(ghostLayer)));
    }
}
