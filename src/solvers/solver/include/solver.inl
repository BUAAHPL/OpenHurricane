#include "solver.hpp"
/*!
 * \file solver.inl
 * \brief The In-Line functions of the <i>solver.hpp</i> file.
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

inline void OpenHurricane::solver::setUsualSolver() noexcept {
    timeType_ = timeTypes::UsualSolver;
}

inline void OpenHurricane::solver::setBDFUnsteadySolver() noexcept {
    timeType_ = timeTypes::BDFUnsteadySolver;
}

inline bool OpenHurricane::solver::hasTPLmtCells() const noexcept {
    return hasTPLmtCells_;
}

inline OpenHurricane::sourceTerms &OpenHurricane::solver::sorcTerm() noexcept {
    return *sorcTermPtr_;
}

inline const OpenHurricane::sourceTerms &OpenHurricane::solver::sorcTerm() const noexcept {
    return *sorcTermPtr_;
}

template <class Type>
inline void OpenHurricane::solver::deletePointer(geometryArray<Type, cellMesh> *_gPtr) const {
    if (_gPtr != nullptr) {
        delete _gPtr;
        _gPtr = nullptr;
    }
}

template <class Type, class GeoMesh>
void OpenHurricane::solver::extrapolate(const Array<Array<Type>> &faceVar,
                                        geometryArray<Type, GeoMesh> &var) {
    const registerTable &tb = iter_.tb();
    const runtimeMesh &mesh = tb.findObject<runtimeMesh>("CFDIter");
    const faceList &fL = mesh.faces();
    const faceZoneList &fZL = mesh.faceZones();
    //const integerListList& SNC = mesh.secondNeighbourCells();

    for (integer fz = 0; fz < faceVar.size(); fz++) {
        for (integer i = 0; i < faceVar[fz].size(); i++) {
            integer index = fZL[fz].firstIndex();
            integer inter = fL[i + index].leftCell();
            integer ghost = fL[i + index].rightCell();
            var[ghost] = real(2.0) * faceVar[fz][i] - var[inter];
        }
    }

    const cutZoneList &cZL = mesh.cutZones();
    integer ghostLayer = mesh.ghostCellsType();
    integer czs = cZL.size() / ghostLayer;

    for (integer i = 0; i < czs; i++) {
        cZL[i].transfer(var);
    }

    for (integer layerI = 1; layerI < mesh.ghostCellsType(); layerI++) {        
        for (integer i = layerI * czs; i < (layerI + 1) * czs; i++) {
            cZL[i].transfer(var);
        }
    }

    /*if (mesh.ghostCellsType() > 1)
    {
            const integerListList& SNC = mesh.secondNeighbourCells();
            for (integer fz = 0; fz < faceVar.size(); fz++)
            {
                    for (integer i = 0; i < faceVar[fz].size(); i++)
                    {
                            for (integer layer = 1; layer < mesh.ghostCellsType(); layer++)
                            {
                                    integer index = fZL[fz].firstIndex();
                                    integer inter = SNC[i + index][2 * layer - 2];
                                    integer ghost = SNC[i + index][2 * layer - 1];
                                    var[ghost] = real(2.0) * faceVar[fz][i] - var[inter];
                            }
                    }
            }
    }*/
}
inline void OpenHurricane::solver::Ac(const integer celli, const vector &normal,
                                      realSquareMatrix &AC) const {
    Ac(celli, normal, real(0), AC);
}
inline OpenHurricane::solver::~solver() noexcept {
    clear();
}

inline void OpenHurricane::solver::solve() {
    switch (timeType_) {
    case OpenHurricane::solver::timeTypes::UsualSolver:
        solving();
        break;
    case OpenHurricane::solver::timeTypes::BDFUnsteadySolver:
        BDFSolve();
        break;
    default:
        solving();
        break;
    }
}

hur_nodiscard inline const OpenHurricane::runtimeMesh &
OpenHurricane::solver::mesh() const noexcept {
    return flowPtr_->mesh();
}

namespace OpenHurricane {
    inline void OpenHurricane::solver::initialize() {
        //const_cast<runtimeMesh&>(mesh()).initResidualsList(iter().residualsNameList());
        //iter_.setMyMonitorPtr(new monitors(iter_, mesh()));
        iter_.initializing(*flowPtr_);
    }

    inline void OpenHurricane::solver::iterRefresh() {       
        flowPtr_->CFLFlag() = 1;
    }

    hur_nodiscard inline const iteration &OpenHurricane::solver::iter() const noexcept {
        return iter_;
    }

    inline void OpenHurricane::solver::write() {
        iter_.write();
    }

    hur_nodiscard inline realArray &solver::dt() noexcept {
        return dt_;
    }
    hur_nodiscard inline spatialScheme &solver::invFlux() noexcept {
        return *invFluxPtr_;
    }
    hur_nodiscard inline timeMarching &solver::marching() noexcept {
        return *timeMarcingPtr_;
    }
    hur_nodiscard inline cellRealArray &solver::shockFactor() noexcept {
        return flowPtr_->shockFactor();
    }
    hur_nodiscard inline const cellRealArray &solver::shockFactor() const noexcept {
        return flowPtr_->shockFactor();
    }
    hur_nodiscard inline cellRealArray &solver::rho() noexcept {
        return flowPtr_->rho();
    }
    hur_nodiscard inline const cellRealArray &solver::rho() const noexcept {
        return flowPtr_->rho();
    }
    hur_nodiscard inline cellRealArray &solver::p() noexcept {
        return flowPtr_->p();
    }
    hur_nodiscard inline const cellRealArray &solver::p() const noexcept {
        return flowPtr_->p();
    }
    hur_nodiscard inline cellRealArray &solver::E() noexcept {
        return flowPtr_->E();
    }
    hur_nodiscard inline const cellRealArray &solver::E() const noexcept {
        return flowPtr_->E();
    }
    hur_nodiscard inline cellRealArray &solver::T() noexcept {
        return flowPtr_->T();
    }
    hur_nodiscard inline const cellRealArray &solver::T() const noexcept {
        return flowPtr_->T();
    }
    hur_nodiscard inline cellRealArray &solver::mu() noexcept {
        return flowPtr_->mul();
    }
    hur_nodiscard inline const cellRealArray &solver::mu() const noexcept {
        return flowPtr_->mul();
    }
    hur_nodiscard inline cellRealArray &solver::mut() noexcept {
        return flowPtr_->mut();
    }
    hur_nodiscard inline const cellRealArray &solver::mut() const noexcept {
        return flowPtr_->mut();
    }
    hur_nodiscard inline mixture &solver::mixtures() noexcept {
        return flowPtr_->mixtures();
    }
    hur_nodiscard inline const mixture &solver::mixtures() const noexcept {
        return flowPtr_->mixtures();
    }

    hur_nodiscard inline speciesList &solver::specTable() noexcept {
        return flowPtr_->mixtures().species();
    }
    hur_nodiscard inline const speciesList &solver::specTable() const noexcept {
        return flowPtr_->mixtures().species();
    }
    hur_nodiscard inline cellVectorArray &solver::v() noexcept {
        return flowPtr_->v();
    }
    hur_nodiscard inline const cellVectorArray &solver::v() const noexcept {
        return flowPtr_->v();
    }

    hur_nodiscard inline OpenHurricane::real OpenHurricane::solver::prl() const noexcept {
        return flowPtr_->Prl();
    }

    hur_nodiscard inline real solver::prt() const noexcept {
        return flowPtr_->Prt();
    }

    hur_nodiscard inline real solver::sct() const noexcept {
        return flowPtr_->Sct();
    }

    hur_nodiscard inline cellRealArray &OpenHurricane::solver::kappal() noexcept {
        return flowPtr_->kappal();
    }

    hur_nodiscard inline const cellRealArray &solver::kappal() const noexcept {
        return flowPtr_->kappal();
    }

    hur_nodiscard inline cellRealArray &solver::kappat() noexcept {
        return flowPtr_->kappat();
    }

    hur_nodiscard inline const cellRealArray &solver::kappat() const noexcept {
        return flowPtr_->kappat();
    }

    hur_nodiscard inline cellRealArray solver::kappaEff() const {
        return flowPtr_->kappaEff();
    }

    hur_nodiscard inline cellRealArray &solver::cp() noexcept {
        return flowPtr_->thermo().cp();
    }

    hur_nodiscard inline const cellRealArray &solver::cp() const noexcept {
        return flowPtr_->thermo().cp();
    }

    hur_nodiscard inline cellRealArray &solver::gama() noexcept {
        return flowPtr_->gama();
    }

    hur_nodiscard inline const cellRealArray &solver::gama() const noexcept {
        return flowPtr_->gama();
    }

    hur_nodiscard inline PtrList<cellRealArray> &solver::yi() noexcept {
        return mixtures().Yi();
    }

    hur_nodiscard inline const PtrList<cellRealArray> &solver::yi() const noexcept {
        return mixtures().Yi();
    }

    hur_nodiscard inline PtrList<cellRealArray> &solver::hi() noexcept {
        return mixtures().hi();
    }

    hur_nodiscard inline const PtrList<cellRealArray> &solver::hi() const noexcept {
        return mixtures().hi();
    }

    hur_nodiscard inline PtrList<cellRealArray> &solver::Diff() noexcept {
        return flowPtr_->mixtures().Dim();
    }

    hur_nodiscard inline const PtrList<cellRealArray> &
    OpenHurricane::solver::Diff() const noexcept {
        return flowPtr_->mixtures().Dim();
    }

    hur_nodiscard inline integerArray &OpenHurricane::solver::temperatureFlag() noexcept {
        return flowPtr_->temperatureFlag();
    }

    hur_nodiscard inline const integerArray &
    OpenHurricane::solver::temperatureFlag() const noexcept {
        return flowPtr_->temperatureFlag();
    }

    hur_nodiscard inline const integerArray &solver::pressureFlag() const noexcept {
        return flowPtr_->pressureFlag();
    }

    hur_nodiscard inline integerArray &solver::pressureFlag() noexcept {
        return flowPtr_->pressureFlag();
    }

    hur_nodiscard inline const cellIntegerArray &solver::CFLFlag() const noexcept {
        return flowPtr_->CFLFlag();
    }
    hur_nodiscard inline cellIntegerArray &solver::CFLFlag() noexcept {
        return flowPtr_->CFLFlag();
    }

    inline void solver::firstSolveSplittingSource(const real phyDt) {
        // Do nothing
        // Only implemented in operator-splitting solver
    }
    inline void solver::secondSolveSplittingSource(const real phyDt) {
        // Do nothing
        // Only implemented in operator-splitting solver
    }
    inline void solver::previousSource() {}
    inline void solver::postSource() {}
   
} // namespace OpenHurricane

template <class Type>
hur_nodiscard OpenHurricane::Array<Type>
OpenHurricane::solver::gatherCellDataInMaster(const Array<Type> &cellQ) const {
    integerList nSizeL;
    integerList displs;
    integer allSize = 0;

    Array<Type> rootF;
    if (HurMPI::parRun()) {
        nSizeL.resize(HurMPI::getProcSize(), Zero);
        nSizeL[HurMPI::getProcRank()] = mesh().nCells();
        if (HurMPI::master()) {
            displs.resize(HurMPI::getProcSize(), Zero);
        }
        HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());
        allSize = 0;
        if (HurMPI::master()) {
            for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
            }
            for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                allSize += nSizeL[ip];
            }
        }

        if (HurMPI::master()) {
            rootF.resize(allSize);
        }

        HurMPI::barrier(HurMPI::getComm());
        /*HurMPI::gatherv
        (
            cellQ.data(),
            mesh().nCells(),
            feature<Type>::MPIType,
            rootF.data(),
            nSizeL.data(),
            displs.data(),
            feature<Type>::MPIType,
            HurMPI::masterNo(),
            HurMPI::getComm()
        );*/

        HurMPI::Request request;
        HurMPI::igatherv(cellQ.data(), mesh().nCells(), feature<Type>::MPIType, rootF.data(),
                         nSizeL.data(), displs.data(), feature<Type>::MPIType, HurMPI::masterNo(),
                         HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);
    } else {
        rootF.resize(mesh().nCells());
        for (integer i = 0; i < mesh().nCells(); i++) {
            rootF[i] = cellQ[i];
        }
    }
    return rootF;
}

hur_nodiscard inline bool OpenHurricane::solver::useLowMachPrecon() const noexcept {
    return useLowMachPrecon_;
}

inline OpenHurricane::cellRealArray &OpenHurricane::solver::a4() {
    if (a4Ptr_ == nullptr) {
        a4Ptr_ = new cellRealArray(object("a4ForLowMachPrecond", mesh(), object::NOT_WRITE), mesh(),
                                   Zero);
    }
    return *a4Ptr_;
}

inline const OpenHurricane::cellRealArray &OpenHurricane::solver::a4() const {
    if (a4Ptr_ == nullptr) {
        a4Ptr_ = new cellRealArray(object("a4ForLowMachPrecond", mesh(), object::NOT_WRITE), mesh(),
                                   Zero);
    }
    return *a4Ptr_;
}

inline OpenHurricane::cellRealArray &OpenHurricane::solver::a5() {
    if (a5Ptr_ == nullptr) {
        a5Ptr_ = new cellRealArray(object("a5ForLowMachPrecond", mesh(), object::NOT_WRITE), mesh(),
                                   Zero);
    }
    return *a5Ptr_;
}

inline const OpenHurricane::cellRealArray &OpenHurricane::solver::a5() const {
    if (a5Ptr_ == nullptr) {
        a5Ptr_ = new cellRealArray(object("a5ForLowMachPrecond", mesh(), object::NOT_WRITE), mesh(),
                                   Zero);
    }
    return *a5Ptr_;
}