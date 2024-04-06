/*!
 * \file reordering.cpp
 * \brief Main subroutines for reordering mesh.
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
#include "reordering.hpp"
#include "HurMPI.hpp"

 namespace OpenHurricane{
	createClassNameStr(reordering,"reordering");
}

namespace OpenHurricane {
     createObjFty(reordering,controller);
}

OpenHurricane::uniquePtr<OpenHurricane::reordering>
OpenHurricane::reordering::creator(decomposingMeshByMetis &dcpM, const controller &cont) {

    string solverType = cont.findWord(reordering::className_);

    
	defineInObjCreator(reordering,static_cast<std::string>(solverType),controller,(dcpM));
}

void OpenHurricane::reordering::computeArray() {
    integer ints = 0;
    for (integer zoneI = 0; zoneI < decomposedMesh_.faceZones().size(); ++zoneI) {
        if (decomposedMesh_.faceZones()[zoneI].isInterior() &&
            (!decomposedMesh_.faceZones()[zoneI].isInterface())) {
            ints += decomposedMesh_.faceZones()[zoneI].size();
        }
    }

    nvtxs_ = decomposedMesh_.cells().size();
    xadj_.resize(nvtxs_ + 1);
    adjncy_.resize(2 * ints);
    perm_.resize(nvtxs_);
    iperm_.resize(nvtxs_);

    /*
     * The adjacency structure of the graph is stored as follows. Assuming that vertex numbering starts from 0 (C style),
     * then the adjacency list of vertex i is stored in array adjncy starting at index xadj[i] and ending at (but not
     * including) index xadj[i+1] (i.e., adjncy[xadj[i]] through and including adjncy[xadj[i+1]-1]). That
     * is, for each vertex i, its adjacency list is stored in consecutive locations in the array adjncy, and the array xadj
     * is used to point to where it begins and where it ends. Figure (b) illustrates the CSR format for the 15-vertex graph
     * shown in Figure (a).
     *								0——1——2——3——4
     *								|  |  |  |  |
     *								5——6——7——8——9
     *								|  |  |  |  |
     *								10—11—12—13—14
     *                            (a) A sample graph
     *  __________________________________________________________________________________________________________________
     * |   xadj     0 2 5 8 11 13 16 20 24 28 31 33 36 39 42 44                                                           |
     * | adjncy     1 5 0 2 6 1 3 7 2 4 8 3 9 0 6 10 1 5 7 11 2 6 8 12 3 7 9 13 4 8 14 5 11 6 10 12 7 11 13 8 12 14 9 13  |
     *  ——————————————————————————————————————————————————————————————————————————————————————————————————————————————————
     *                            (b) CSR format
     */
    integer k = 0;
    xadj_[0] = 0;
    for (integer celli = 0; celli < nvtxs_; celli++) {
        xadj_[celli + 1] = xadj_[celli];

        for (integer fci = 0; fci < decomposedMesh_.cells()[celli].faceSize(); fci++) {
            integer facei = decomposedMesh_.cells()[celli].facei(fci);
            const face &fc = decomposedMesh_.faces()[facei];

            if (fc.isInterior() && (!fc.isInterface())) {
                adjncy_[k] = fc.leftCell() + fc.rightCell() - celli;
                xadj_[celli + 1] += 1;
                k += 1;
            }
        }
    }
}

void OpenHurricane::reordering::update() {
    cellList cellTemp(nvtxs_);
    for (integer n = 0; n < nvtxs_; n++) {
        integer i = iperm_[n];
        perm_[i] = n;
        //cellTemp[n] = decomposedMesh_.cells()[i];
        cellTemp[n].transfer(decomposedMesh_.cells()[i]);
    }

    /*decomposedMesh_.cells() = cellTemp;
    cellTemp.clear();*/

    decomposedMesh_.cells().transfer(cellTemp);

    /*!\brief Update cell order for face structure.*/
    for (integer i = 0; i < decomposedMesh_.faces().size(); ++i) {
        integer n = decomposedMesh_.faces()[i].leftCell();
        decomposedMesh_.faces()[i].leftCell() = perm_[n];
        n = decomposedMesh_.faces()[i].rightCell();
        if (n != -1) {
            decomposedMesh_.faces()[i].rightCell() = perm_[n];
        }
    }

    updateFacePairMap();

    updatePerPairMap();

    /*!\brief Update cell order for interface structure(undone).*/
}

void OpenHurricane::reordering::updateFacePairMap() {
    if (decomposedMesh_.facePairMaps().size() != 0) {
        std::map<integer, integerArray> &facePairMaps = decomposedMesh_.facePairMaps();
        const faceList &fL = decomposedMesh_.faces();
        integerArray size(HurMPI::getProcSize(), Zero);
        for (std::map<integer, integerArray>::const_iterator iter = facePairMaps.begin();
             iter != facePairMaps.end(); ++iter) {
            integer ip = iter->second[0];
            size[ip]++;
        }

        integerArrayArray faceSendBuf(HurMPI::getProcSize());
        integerArrayArray cellSendBuf(HurMPI::getProcSize());
        integerArrayArray faceRecvBuf(HurMPI::getProcSize());
        integerArrayArray cellRecvBuf(HurMPI::getProcSize());
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            faceSendBuf[ip].resize(size[ip]);
            cellSendBuf[ip].resize(size[ip]);
            faceRecvBuf[ip].resize(size[ip]);
            cellRecvBuf[ip].resize(size[ip]);
            size[ip] = Zero;
        }

        for (std::map<integer, integerArray>::const_iterator iter = facePairMaps.begin();
             iter != facePairMaps.end(); ++iter) {
            integer ip = iter->second[0];
            integer faceI = iter->second[2];
            integer cellI = fL[iter->first].leftCell();
            faceSendBuf[ip][size[ip]] = faceI;
            cellSendBuf[ip][size[ip]] = cellI;
            size[ip]++;
        }

        integer countIP = Zero;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            if (size[ip] != 0) {
                countIP++;
            }
        }

        List<HurMPI::Request> request(countIP * 2);
        integer countReq = 0;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            if (size[ip] != 0) {
                HurMPI::irecv(faceRecvBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 0,
                              HurMPI::getComm(), &request[countReq]);
                countReq++;
            }
        }
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            if (size[ip] != 0) {
                HurMPI::isend(faceSendBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 0,
                              HurMPI::getComm(), &request[countReq]);
                countReq++;
            }
        }
        HurMPI::waitall(request.size(), request.data(), MPI_STATUSES_IGNORE);
        countReq = 0;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            if (size[ip] != 0 && !HurMPI::isThisProc(ip)) {
                HurMPI::irecv(cellRecvBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 1,
                              HurMPI::getComm(), &request[countReq]);
                countReq++;
            }
        }
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            if (size[ip] != 0 && !HurMPI::isThisProc(ip)) {
                HurMPI::isend(cellSendBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 1,
                              HurMPI::getComm(), &request[countReq]);
                countReq++;
            }
        }
        HurMPI::waitall(request.size(), request.data(), MPI_STATUSES_IGNORE);

        // 下面用法在某些情况下会造成程序死锁
        /*for (integer ip = 0; ip < HurMPI::getProcSize(); ip++)
        {
                if (size[ip] != 0)
                {
                        HurMPI::send(faceSendBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 0, HurMPI::getComm());
                        HurMPI::send(cellSendBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 1, HurMPI::getComm());
                        HurMPI::recv(faceRecvBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 0, HurMPI::getComm(), &status);
                        HurMPI::recv(cellRecvBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 1, HurMPI::getComm(), &status);
                }
        }*/

        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            for (integer i = 0; i < size[ip]; i++) {
                integer faceI = faceRecvBuf[ip][i];
                facePairMaps[faceI][3] = cellRecvBuf[ip][i];
            }
        }
    }
}

void OpenHurricane::reordering::updatePerPairMap() {
    if (decomposedMesh_.perPairMaps().size() != 0) {
        std::map<integer, integerArray> &perPairMaps = decomposedMesh_.perPairMaps();
        const faceList &fL = decomposedMesh_.faces();
        integerArray size(HurMPI::getProcSize(), Zero);

        integerArrayArray faceSendBuf(HurMPI::getProcSize());
        integerArrayArray cellSendBuf(HurMPI::getProcSize());
        integerArrayArray faceRecvBuf(HurMPI::getProcSize());
        integerArrayArray cellRecvBuf(HurMPI::getProcSize());
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            faceSendBuf[ip].resize(static_cast<integer>(perPairMaps.size()));
            cellSendBuf[ip].resize(static_cast<integer>(perPairMaps.size()));
            faceRecvBuf[ip].resize(static_cast<integer>(perPairMaps.size()));
            cellRecvBuf[ip].resize(static_cast<integer>(perPairMaps.size()));
            size[ip] = 0;
        }

        for (std::map<integer, integerArray>::const_iterator iter = perPairMaps.begin();
             iter != perPairMaps.end(); ++iter) {
            integer ip = iter->second[0];
            integer faceI = iter->second[1];
            integer cellI = fL[iter->first].leftCell();
            if (HurMPI::isThisProc(ip)) {
                perPairMaps[faceI][2] = cellI;
            } else {
                faceSendBuf[ip][size[ip]] = faceI;
                cellSendBuf[ip][size[ip]] = cellI;
                size[ip]++;
            }
        }
        integer countIP = Zero;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            if (size[ip] != 0) {
                countIP++;
            }
        }

        List<HurMPI::Request> request(countIP * 2);

        integer countReq = 0;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            if (size[ip] != 0) {
                HurMPI::irecv(faceRecvBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 0,
                              HurMPI::getComm(), &request[countReq]);
                countReq++;
            }
        }

        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            if (size[ip] != 0) {
                HurMPI::isend(faceSendBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 0,
                              HurMPI::getComm(), &request[countReq]);
                countReq++;
            }
        }

        HurMPI::waitall(countIP * 2, request.data(), MPI_STATUSES_IGNORE);

        countReq = 0;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            if (size[ip] != 0) {
                HurMPI::irecv(cellRecvBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 1,
                              HurMPI::getComm(), &request[countReq]);
                countReq++;
            }
        }

        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            if (size[ip] != 0) {
                HurMPI::isend(cellSendBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 1,
                              HurMPI::getComm(), &request[countReq]);
                countReq++;
            }
        }
        HurMPI::waitall(countIP * 2, request.data(), MPI_STATUSES_IGNORE);

        // 下面用法在某些情况下会造成程序死锁
        /*HurMPI::Status status;
        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++)
        {
                if (size[ip] != 0)
                {
                        HurMPI::send(faceSendBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 0, HurMPI::getComm());
                        HurMPI::send(cellSendBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 1, HurMPI::getComm());
                        HurMPI::recv(faceRecvBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 0, HurMPI::getComm(), &status);
                        HurMPI::recv(cellRecvBuf[ip].data(), size[ip], feature<integer>::MPIType, ip, 1, HurMPI::getComm(), &status);
                }
        }*/

        for (integer ip = 0; ip < HurMPI::getProcSize(); ip++) {
            for (integer i = 0; i < size[ip]; i++) {
                integer faceI = faceRecvBuf[ip][i];
                perPairMaps[faceI][2] = cellRecvBuf[ip][i];
            }
        }
    }
}

OpenHurricane::integer OpenHurricane::reordering::getBandwidth() const {
    integer nb = 0;
    for (integer n = 0; n < decomposedMesh_.cells().size(); ++n) {
        integer cmax = n;
        integer cmin = n;
        for (integer fli = 0; fli < decomposedMesh_.cells()[n].faceSize(); ++fli) {
            const integer fi = decomposedMesh_.cells()[n].facei(fli);
            integer k = decomposedMesh_.faces()[fi].leftCell() +
                        decomposedMesh_.faces()[fi].rightCell() - n;
            if (k >= 0) {
                cmax = max(cmax, k);
                cmin = min(cmin, k);
            }
        }
        nb = max(nb, (cmax - cmin + 1));
    }
    return nb;
}

void OpenHurricane::reordering::printBandwidthReduction(const integer bw1, const integer bw2) const {
    /*integer bnndw1 = bw1;
    integer bnndw2 = bw2;
    HurMPI::reduce(bnndw1, MPI_SUM);
    HurMPI::reduce(bnndw2, MPI_SUM);
    Pout << "      bandwidth reduced to " << bnndw1 << "/" << bnndw2 << " = "
            << real(bnndw1) / real(bnndw2) << std::endl;*/
    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        if (HurMPI::isThisProc(ip)) {
            std::cout << "      bandwidth reduced to " << bw1 << "/" << bw2 << " = "
                      << real(bw1) / real(bw2) << " in process: " << ip << std::endl;
        }
        HurMPI::barrier();
    }
}
