/*!
 * \file searchProcedures.cpp
 * \brief Main subroutines for search procedures of distance method.
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

#include "searchProcedures.hpp"
#include "HurMPI.hpp"
#include "Lists.hpp"
OpenHurricane::faceZonePackage::faceZonePackage() : faceCentre_(), faceArea_(), fp_() {}

namespace OpenHurricane {
    createClassNameStr(searchProcedures, "searchProcedures");
    registerObjFty(distanceMethod, searchProcedures, controller);
} // namespace OpenHurricane

void OpenHurricane::searchProcedures::setFacePack() {
    const faceZoneList &fzl = mesh_.faceZones();
    const auto &fA = mesh_.faceArea();
    const auto &fC = mesh_.faceCentre();

    integerList faceSize(HurMPI::getProcSize(), Zero);
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const integer fzi = zoneIdList_[i];
        faceSize[HurMPI::getProcRank()] += fzl[fzi].size();
    }

    integer fid = 0;
    faceIdList_.resize(faceSize[HurMPI::getProcRank()]);
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const integer fzi = zoneIdList_[i];
        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            faceIdList_[fid++] = fi;
        }
    }

    HurMPI::barrier(HurMPI::getComm());
    HurMPI::allGatherList(faceSize, HurMPI::getComm());
    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        facePack_[ip].faceArea().resize(faceSize[ip]);
        facePack_[ip].faceCentre().resize(faceSize[ip]);
        facePack_[ip].fp().resize(faceSize[ip]);
    }

    for (integer i = 0; i < faceSize.size(); ++i) {
        nFaces_ += faceSize[i];
    }

    integer fsize = 0;
    for (integer i = 0; i < zoneIdList_.size(); ++i) {
        const integer fzi = zoneIdList_[i];
        const integer ip = HurMPI::getProcRank();
        for (integer fi = fzl[fzi].firstIndex(); fi <= fzl[fzi].lastIndex(); ++fi) {
            facePack_[ip].faceArea()[fsize] = fA[fi];
            facePack_[ip].faceCentre()[fsize++] = fC[fi];
        }
    }

    for (int ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        integer num = 6 * faceSize[ip];
        real *send = new real[num];
        if (HurMPI::isThisProc(ip)) {
            for (integer i = 0; i < faceSize[ip]; i++) {
                send[i * 6] = facePack_[ip].faceArea()[i][0];
                send[i * 6 + 1] = facePack_[ip].faceArea()[i][1];
                send[i * 6 + 2] = facePack_[ip].faceArea()[i][2];
                send[i * 6 + 3] = facePack_[ip].faceCentre()[i][0];
                send[i * 6 + 4] = facePack_[ip].faceCentre()[i][1];
                send[i * 6 + 5] = facePack_[ip].faceCentre()[i][2];
            }
        }
        HurMPI::bcast(send, num, feature<real>::MPIType, ip, HurMPI::getComm());
        if (!HurMPI::isThisProc(ip)) {
            for (integer i = 0; i < faceSize[ip]; i++) {
                facePack_[ip].faceArea()[i][0] = send[i * 6];
                facePack_[ip].faceArea()[i][1] = send[i * 6 + 1];
                facePack_[ip].faceArea()[i][2] = send[i * 6 + 2];
                facePack_[ip].faceCentre()[i][0] = send[i * 6 + 3];
                facePack_[ip].faceCentre()[i][1] = send[i * 6 + 4];
                facePack_[ip].faceCentre()[i][2] = send[i * 6 + 5];
            }
        }
        HurDeleteDynArray(send);
        HurMPI::barrier(HurMPI::getComm());
    }
    setFacePoint();
}

void OpenHurricane::searchProcedures::setFacePoint() {
    const faceZoneList &fzl = mesh_.faceZones();
    const auto &fA = mesh_.faceArea();
    const auto &fC = mesh_.faceCentre();
    const auto &cS = mesh_.cells();
    const auto &fS = mesh_.faces();
    const auto &ps = mesh_.points();

    integerList adj(faceIdList_.size() + 1, Zero);
    for (integer i = 0; i < faceIdList_.size(); ++i) {
        adj[i + 1] = adj[i] + fS[faceIdList_[i]].size();
    }

    integerList nSizeL(HurMPI::getProcSize(), Zero);
    integerList displs;
    displs.resize(HurMPI::getProcSize(), Zero);

    nSizeL[HurMPI::getProcRank()] = faceIdList_.size();
    HurMPI::allGatherList(nSizeL, HurMPI::getComm());

    for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
        displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
    }

    integer ssi = 0;
    for (integer i = 0; i < nSizeL.size(); ++i) {
        ssi += nSizeL[i];
    }

    integerList allAdj(ssi);
    HurMPI::Request request;
    HurMPI::iallGatherv(adj.data(), adj.size() - 1, feature<integer>::MPIType, allAdj.data(),
                        nSizeL.data(), displs.data(), feature<integer>::MPIType, HurMPI::getComm(),
                        &request);
    HurMPI::wait(&request, MPI_STATUSES_IGNORE);

    integerList nSizeOfP(HurMPI::getProcSize(), Zero);
    nSizeOfP[HurMPI::getProcRank()] = adj.last();
    HurMPI::allGatherList(nSizeOfP, HurMPI::getComm());
    displs[0] = 0;
    for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
        displs[ip] = displs[ip - 1] + nSizeOfP[ip - 1];
    }

    realArray px(adj.last());
    realArray py(adj.last());
    realArray pz(adj.last());
    integer count = Zero;
    for (integer i = 0; i < faceIdList_.size(); ++i) {
        for (integer pj = 0; pj < fS[faceIdList_[i]].size(); ++pj) {
            px[count] = ps[fS[faceIdList_[i]][pj]].x();
            py[count] = ps[fS[faceIdList_[i]][pj]].y();
            pz[count] = ps[fS[faceIdList_[i]][pj]].z();
            count++;
        }
    }
    integer ssxyz = 0;
    for (integer i = 0; i < nSizeL.size(); ++i) {
        ssxyz += nSizeOfP[i];
    }
    realArray allPx(ssxyz);
    realArray allPy(ssxyz);
    realArray allPz(ssxyz);

    List<HurMPI::Request> requestxyz(3);
    HurMPI::iallGatherv(px.data(), px.size(), feature<real>::MPIType, allPx.data(), nSizeOfP.data(),
                        displs.data(), feature<real>::MPIType, HurMPI::getComm(), &requestxyz[0]);

    HurMPI::iallGatherv(py.data(), py.size(), feature<real>::MPIType, allPy.data(), nSizeOfP.data(),
                        displs.data(), feature<real>::MPIType, HurMPI::getComm(), &requestxyz[1]);

    HurMPI::iallGatherv(pz.data(), pz.size(), feature<real>::MPIType, allPz.data(), nSizeOfP.data(),
                        displs.data(), feature<real>::MPIType, HurMPI::getComm(), &requestxyz[2]);
    HurMPI::waitall(requestxyz.size(), requestxyz.data(), MPI_STATUSES_IGNORE);

    integer offsetf = Zero;
    integer offsetp = Zero;
    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        for (integer i = 0; i < facePack_[ip].faceArea().size(); ++i) {
            const integer begin = allAdj[i + offsetf];
            integer end;
            if (i == facePack_[ip].faceArea().size() - 1) {
                end = nSizeOfP[ip] + offsetp;
            } else {
                end = allAdj[i + 1 + offsetf];
            }

            integer counti = 0;
            facePack_[ip].fp()[i].resize(end - begin);
            for (integer pi = begin; pi < end; ++pi) {
                facePack_[ip].fp()[i][counti].x() = allPx[pi];
                facePack_[ip].fp()[i][counti].y() = allPy[pi];
                facePack_[ip].fp()[i][counti].z() = allPz[pi];
                counti++;
            }
        }
        offsetp += nSizeOfP[ip];
        offsetf += nSizeL[ip];
    }
}

void OpenHurricane::searchProcedures::distTransfer() {
    const auto &cutZoneL = mesh_.cutZones();
    for (integer cuti = 0; cuti < cutZoneL.size(); ++cuti) {
        cutZoneL[cuti].transfer(iw1_);
        cutZoneL[cuti].transfer(iw2_);
    }
}

void OpenHurricane::searchProcedures::smooth(cellRealArray &y, integer &t) {
    t = 0;
    const integer nCell = mesh_.nCells();
    const auto &cC = mesh_.cellCentre();
    const auto &cS = mesh_.cells();
    const auto &fS = mesh_.faces();

    for (integer n = 0; n < nCell; ++n) {
        iw1_[n] = iw2_[n];
        for (integer i = 0; i < cS[n].faceSize(); ++i) {
            const integer fi = cS[n].facei(i);
            integer m = fS[fi].leftCell() + fS[fi].rightCell() - n;
            integer p = iw2_[m];
            if (p != -1) {
                integer ip = link_[p][0];
                integer ii = link_[p][1];
                vector cw = facePack_[ip].faceCentre()[ii] - cC[n];
                real r2 = cw.magnitude();
                if (r2 < y[n]) {
                    y[n] = r2;
                    iw2_[n] = p;
                }
            }
        }
        if (iw2_[n] != iw1_[n]) {
            t++;
        }
    }
}

OpenHurricane::real OpenHurricane::searchProcedures::calcDist(const integer ip, const integer ii,
                                                      const vector &cC) const {
    const vector &nn = facePack_[ip].faceArea()[ii];
    const vector &rr = facePack_[ip].faceCentre()[ii];
    real r = 0;
    const auto pInFace = faceCalc::projectionInface(cC, rr, nn);
    bool isNormal = false;
    if (facePack_[ip].fp()[ii].size() == 3) {
        isNormal =
            faceCalc::pointInTriangularFace(pInFace, facePack_[ip].fp()[ii][0],
                                            facePack_[ip].fp()[ii][1], facePack_[ip].fp()[ii][2]);
    } else if (facePack_[ip].fp()[ii].size() == 4) {
        isNormal = faceCalc::pointInQuadrilateralFace(
            pInFace, facePack_[ip].fp()[ii][0], facePack_[ip].fp()[ii][1],
            facePack_[ip].fp()[ii][2], facePack_[ip].fp()[ii][3]);
    }

    if (isNormal) {
        r = fabs((rr - cC) * nn.normalized());
    } else {
        r = (rr - cC).magnitude();
    }
    return r;
}

OpenHurricane::searchProcedures::searchProcedures(const controller &cont, const runtimeMesh &mesh,
                                              const integerList zoneIdL)
    : distanceMethod(cont, mesh, zoneIdL), facePack_(HurMPI::getProcSize()), faceIdList_(), iw1_(),
      iw2_(), link_(), itt_(), nFaces_(0), normal_(false) {
    setFacePack();
}

OpenHurricane::searchProcedures::searchProcedures(const runtimeMesh &mesh, const integerList zoneIdL)
    : distanceMethod(mesh, zoneIdL), facePack_(HurMPI::getProcSize()), faceIdList_(), iw1_(),
      iw2_(), link_(), itt_(), nFaces_(0), normal_(false) {
    setFacePack();
}

bool OpenHurricane::searchProcedures::getDistance(cellRealArray &y, cellVectorArray &ny) {
    const integer nTotalCells = mesh_.nTotalCells();
    const faceList &fs = mesh_.faces();
    const auto &cC = mesh_.cellCentre();
    const integer nCell = mesh_.nCells();

    iw1_.resize(nTotalCells);
    iw2_.resize(nTotalCells);
    link_.resize(nFaces_);
    itt_.resize(HurMPI::getProcSize());

    iw1_ = -1;
    iw2_ = -1;

    y = 1e20;

    integer ipw = 0;
    integer iw = 0;
    integer pw = 0;
    integer p = 0;
    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
        for (integer i = 0; i < facePack_[ip].faceArea().size(); ++i) {
            link_[p][0] = ip;
            link_[p][1] = i;
            if (ipw == 0) {
                ipw = ip;
                iw = i;
                pw = p;
            }
            if (HurMPI::isThisProc(ip)) {
                const auto &cl = fs[faceIdList_[i]].leftCell();
                vector cw = facePack_[ip].faceCentre()[i] - cC[cl];
                real r = cw.magnitude();
                if (r < y[cl]) {
                    y[cl] = r;
                    iw1_[cl] = p;
                    iw2_[cl] = p;
                }
            }
            p++;
        }
    }
    integer t;
    integer mm = 1;
    Pout << "      Begin iteration cycle for search procedures" << std::endl;
    while (mm > 0) {
        integer t_sum = 0;
        distTransfer();
        smooth(y, t);
        HurMPI::barrier(HurMPI::getComm());
        HurMPI::gather(&t, 1, feature<integer>::MPIType, itt_.data(), 1, feature<integer>::MPIType,
                       HurMPI::masterNo(), HurMPI::getComm());
        if (HurMPI::master()) {
            for (integer i = 0; i < itt_.size(); ++i) {
                t_sum += itt_[i];
            }
        }
        HurMPI::bcast(&t_sum, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (mm % 10 == 0) {
            Pout << "      Iteration cycle: " << mm << std::endl;
        }
        if (t_sum == 0) {
            Pout << "      End iteration cycle at: " << mm << std::endl;
            break;
        }
        mm++;
    }

    t = 0;
    for (integer n = 0; n < nCell; ++n) {
        p = iw2_[n];
        if (p == -1) {
            iw2_[n] = pw;
            link_[pw][0] = ipw;
            link_[pw][1] = iw;
            y[n] = 1e30;
            continue;
        }
        integer ip = link_[p][0];
        integer ii = link_[p][1];
        vector &nn = facePack_[ip].faceArea()[ii];
        vector &rr = facePack_[ip].faceCentre()[ii];
        real r = 0;
        //r = calcDist(ip, ii, cC[n]);
        r = fabs((rr - cC[n]) * nn.normalized());
        /*if (normal_)
        {
            r = fabs((rr - cC[n]) * nn.normalized());
        }
        else
        {
            r = (rr - cC[n]).magnitude();
        }*/
        //real r = fabs((rr - cC[n]) * nn.normalized());
        //r = (rr - cC[n]).magnitude();

        if (r >= veryTiny) {
            y[n] = r;
        }
    }
    return true;
}