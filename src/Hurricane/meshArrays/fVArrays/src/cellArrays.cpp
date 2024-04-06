/*!
 * \file cellArrays.cpp
 * \brief Main subroutines for cellArrays.
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
#include "cellArrays.hpp"
#include "boundaries.hpp"
#include "cellMesh.hpp"
#include "hdf5I.hpp"
#include "hdf5O.hpp"
#include "meshInterpolation.hpp"
#include "real.hpp"
#include "relayScatter.hpp"
#include "transfer.hpp"

template <>
void OpenHurricane::geometryArray<OpenHurricane::integer, OpenHurricane::cellMesh>::calcTimeSumPtr(
    const real &dt) const {
    if (!timeSumPtr_) {
        string rn = object::name() + "TimeSum";
        timeSumPtr_.reset(new geometryArray<integer, cellMesh>(
            object(rn, object::tb(), object::NOT_WRITE), mesh_, Zero));
    }
    (*timeSumPtr_) += dt * (*this);
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>::getAveAndMaxRHS() const {
    if (!rhsAvePtr_) {
        rhsAvePtr_.reset(new real());
    }
    if (!rhsMaxPtr_) {
        rhsMaxPtr_.reset(new real());
    }

    real rhsAve = 0.0;
    real zRhs = 0.0;
    real rhsMax = 0.0;

    for (integer i = 0; i < mesh().nCells(); ++i) {
        zRhs = fabs(rhs()[i]) / mesh().cellVol()[i];
        rhsAve += sqr(zRhs);
        rhsMax = max(rhsMax, zRhs);
    }
    realArray rhsAveL(HurMPI::getProcSize(), Zero);
    HurMPI::gather(&rhsAve, 1, feature<real>::MPIType, rhsAveL.data(), 1, feature<real>::MPIType,
                   HurMPI::masterNo(), HurMPI::getComm());

    if (HurMPI::master()) {
        rhsAve = 0.0;
        for (integer i = 0; i < rhsAveL.size(); i++) {
            rhsAve += rhsAveL[i];
        }
        rhsAve /= real(mesh().allMeshCellNumber());
        *rhsAvePtr_ = rhsAve;
    }

    HurMPI::bcast(rhsAvePtr_.get(), 1, feature<real>::MPIType, HurMPI::masterNo());
    HurMPI::allReduce(rhsMax, MPI_MAX);
    *rhsMaxPtr_ = rhsMax;
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>::updateBoundary(
    integer layerI) {
    realTransfer mytransfer(this->mesh(), *this, false, false);

    mytransfer.transferInit(layerI);

    if (!boundariesPtr_) {
        const controller &cont = tb().Iteration().cont();
        boundariesPtr_ = new boundaryCondition<real, cellMesh>(mesh_.faceZones(), *this, cont);
    }
    boundariesPtr_->calculate();
    mytransfer.transferring(layerI);
}

template <>
OpenHurricane::real OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>::initialize() {
    const controller &cont = Iteration().cont();
    if (!cont.found("initialization")) {
        LFatal("Cannot find initialization controller in: %s", cont.name().c_str());
    }
    const auto &initCont = cont.subController("initialization");
    if (!boundariesPtr_) {
        boundariesPtr_ = new boundaryCondition<real, cellMesh>(mesh_.faceZones(), *this, cont);
    }
    if (initCont.found("initFromValue")) {
        const auto &initValueCont = initCont.subController("initFromValue");

        if (initValueCont.found(this->name())) {
            const real initValue = initValueCont.findType<real>(this->name(), real(0.0));
            *this = initValue;
            return initValue;
        }
    }
    for (integer i = 0; i < boundariesPtr_->size(); i++) {
        if ((*boundariesPtr_)[i].boundaryfZ().name() == initCont.findWord("initFromBoundary")) {
            *this = ((*boundariesPtr_)[i].initValue());
            return (*boundariesPtr_)[i].initValue();
        }
    }
    checkWarningStr(("The variable: " + this->name() + " is not initialized."));
    return Zero;
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>::writeOutput(
    fileOsstream &fos) const {
    if (argParse::checkWriteSolution()) {
        bool isok = true;

        for (integer n = 0; n < this->internalArraySize(); ++n) {
            if (isnan(this->operator[](n)) || isinf(this->operator[](n))) {
                isok = false;
                break;
            }
        }
        if (!isok) {
            LFatal("Not a number or infinity of %s", this->name().c_str());
        }
    }

    Array<real>::writeToStream(fos, this->internalArraySize());
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>::writeOutputByMaster(
    fileOsstream &fos) const {
    if (HurMPI::parRun()) {
        const auto outSize = this->internalArraySize();
        integerList nSizeL;
        integerList displs;
        integer allSize = 0;
        nSizeL.resize(HurMPI::getProcSize(), Zero);
        nSizeL[HurMPI::getProcRank()] = outSize;
        if (HurMPI::master()) {
            displs.resize(HurMPI::getProcSize(), Zero);
        }
        HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());

        if (HurMPI::master()) {
            for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
            }
            for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                allSize += nSizeL[ip];
            }
        }
        realArray rootF;
        if (HurMPI::master()) {
            rootF.resize(allSize);
        }
        if (argParse::checkWriteSolution()) {
            bool isok = true;

            for (integer n = 0; n < this->internalArraySize(); ++n) {
                if (isnan(this->operator[](n)) || isinf(this->operator[](n))) {
                    isok = false;
                    break;
                }
            }
            if (!isok) {
                LFatal("Not a number or infinity of %s", this->name().c_str());
            }
        }

        HurMPI::barrier(HurMPI::getComm());
        HurMPI::Request request;
        HurMPI::igatherv(data(), this->internalArraySize(), feature<real>::MPIType, rootF.data(),
                         nSizeL.data(), displs.data(), feature<real>::MPIType, HurMPI::masterNo(),
                         HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);

        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&rootF[0]), rootF.byteSize());
        }
    } else {
        writeOutput(fos);
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>::writeOutput(
    fileOsstream &fos, const integer fzid) const {
    const auto &displs = mesh().globalFaceZoneInfo(fzid).faceDispls();
    const auto &recvcnt = mesh().globalFaceZoneInfo(fzid).faceRecvcnt();
    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();
    const auto fsize = mesh().globalFaceZoneInfo(fzid).fZ().size();
    const auto nTotalFaces = mesh().globalFaceZoneInfo(fzid).totalFaces();

    realArray faceValue(fsize, Zero);

    integer count = 0;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        const integer cl = faces[fi].leftCell();
        const integer cr = faces[fi].rightCell();
        faceValue[count++] = 0.5 * this->operator[](cl) + 0.5 * this->operator[](cr);
    }
    if (!HurMPI::parRun()) {
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            integer i = 0;
            while (i < fsize) {
                sstr << std::setprecision(feature<real>::precision) << faceValue[i++] << " ";
                if (i < fsize) {
                    sstr << std::setprecision(feature<real>::precision) << faceValue[i++] << " ";
                }
                if (i < fsize) {
                    sstr << std::setprecision(feature<real>::precision) << faceValue[i++] << " ";
                }
                sstr << "\n";
            }
            fos.os() << sstr.str().c_str();
        } else {
            if (fsize) {
                fos.write(reinterpret_cast<const char *>(&faceValue[0]), fsize * sizeof(real));
            }
        }
        return;
    }
    realArray writeField;
    if (HurMPI::master()) {
        writeField.resize(nTotalFaces, Zero);
    }

    HurMPI::Request request;
    HurMPI::igatherv(faceValue.data(), fsize, feature<real>::MPIType, writeField.data(),
                     recvcnt.data(), displs.data(), feature<real>::MPIType, HurMPI::masterNo(),
                     HurMPI::getComm(), &request);
    HurMPI::wait(&request, MPI_STATUSES_IGNORE);

    if (HurMPI::master()) {
        if (fos.format() == IOsstream::ASCII_FORMAT) {
            std::stringstream sstr;
            integer i = 0;
            while (i < writeField.size()) {
                sstr << std::setprecision(feature<real>::precision) << writeField[i++] << " ";
                if (i < writeField.size()) {
                    sstr << std::setprecision(feature<real>::precision) << writeField[i++] << " ";
                }
                if (i < writeField.size()) {
                    sstr << std::setprecision(feature<real>::precision) << writeField[i++] << " ";
                }
                sstr << "\n";
            }
            fos.os() << sstr.str().c_str();
        } else {
            if (writeField.size()) {
                fos.write(reinterpret_cast<const char *>(&writeField[0]),
                          writeField.size() * sizeof(real));
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>::writeMinMaxOutput(
    fileOsstream &fos, const integer fzid) const {

    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();

    real minX = veryLarge;
    real maxX = -veryLarge;
    for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
        const integer cl = faces[fi].leftCell();
        const integer cr = faces[fi].rightCell();
        real vv = 0.5 * this->operator[](cl) + 0.5 * this->operator[](cr);
        minX = min(minX, vv);
        maxX = max(maxX, vv);
    }
    if (HurMPI::parRun()) {
        HurMPI::reduce(minX, MPI_MIN);
        HurMPI::reduce(maxX, MPI_MAX);
    }
    HurMPI::barrier(HurMPI::getComm());
    if (HurMPI::master()) {
        fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
        fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>::writeRelay(
    hdf5O &fos, const bool writeLast, const bool writeToGroup) const {
    integerList nSizeL;
    integerList displs;
    integer allSize = 0;
    if (HurMPI::parRun()) {
        nSizeL.resize(HurMPI::getProcSize(), Zero);
        nSizeL[HurMPI::getProcRank()] = this->internalArraySize();
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
        realArray rootF;
        if (HurMPI::master()) {
            rootF.resize(allSize);
        }

        HurMPI::barrier(HurMPI::getComm());

        HurMPI::Request request;
        HurMPI::igatherv(data(), this->internalArraySize(), feature<real>::MPIType, rootF.data(),
                         nSizeL.data(), displs.data(), feature<real>::MPIType, HurMPI::masterNo(),
                         HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);

        if (HurMPI::master()) {
            if (writeLast || writeToGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.write(rootF, rootF.size(), groupName, this->name());
                fos.writeRealAttributeToDataset(rhs0_, "rhs0", groupName, this->name());
            } else {
                fos.write(rootF, rootF.size(), this->name());
                fos.writeRealAttributeToDataset(rhs0_, "rhs0", this->name());
            }
        }
    } else {
        if (writeLast || writeToGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.write(*this, this->internalArraySize(), groupName, this->name());
            fos.writeRealAttributeToDataset(rhs0_, "rhs0", groupName, this->name());
        } else {
            fos.write(*this, this->internalArraySize(), this->name());
            fos.writeRealAttributeToDataset(rhs0_, "rhs0", this->name());
        }
    }
    if (writeLast) {

        if (lastArrayArray_.size() == 0) {
            return;
        }
        if (HurMPI::parRun()) {
            for (integer k = 0; k < lastArrayArray_.size(); ++k) {
                realArray rootF;
                if (HurMPI::master()) {
                    rootF.resize(allSize);
                }
                HurMPI::barrier(HurMPI::getComm());

                HurMPI::Request request;
                HurMPI::igatherv(lastArrayArray_[k].data(), this->internalArraySize(),
                                 feature<real>::MPIType, rootF.data(), nSizeL.data(), displs.data(),
                                 feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm(),
                                 &request);
                HurMPI::wait(&request, MPI_STATUSES_IGNORE);

                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(k + 1);
                    fos.write(rootF, rootF.size(), groupName, this->name());
                }
            }
        } else {
            for (integer k = 0; k < lastArrayArray_.size(); ++k) {
                string groupName = "timeGroupm";
                groupName += toString(k + 1);
                fos.write(lastArrayArray_[k], this->internalArraySize(), groupName, this->name());
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>::readRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {
    if (HurMPI::parRun()) {
        realArray allTmpData;
        realArray relayData;
        if (HurMPI::master()) {
            if (readLast || readFromGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.read(allTmpData, groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_, "rhs0", groupName, this->name());
            } else {
                fos.read(allTmpData, this->name());
                fos.readRealAttributeFromDataset(rhs0_, "rhs0", this->name());
            }
            const std::map<integer, integer> &indexMap = mesh().indexMap();
            const auto &originCellIndex = mesh().originCellIndex();
            relayData.resize(allTmpData.size());

            integer offset = 0;

            for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                    std::map<integer, integer>::const_iterator iter =
                        indexMap.find(originCellIndex[czi][i]);
                    relayData[i + offset] = allTmpData[iter->second];
                }
                offset += originCellIndex[czi].size();
            }
        }

        realArray tmpData;
        relayScatterFunc::relayScatter<real>(mesh().cellOfProc(), relayData, tmpData);
        relayScatterFunc::relayReorder<real>(tmpData, *this, mesh().perm(), mesh().iperm());

        integer ng = 0;
        if (HurMPI::master()) {
            if (!fos.exist("nTimeGroups")) {
                fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
            }
        }
        HurMPI::bcast(&ng, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(i + 1);
                    fos.read(allTmpData, groupName, this->name());
                    const std::map<integer, integer> &indexMap = mesh().indexMap();
                    const auto &originCellIndex = mesh().originCellIndex();

                    integer offset = 0;
                    for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                        for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                            std::map<integer, integer>::const_iterator iter =
                                indexMap.find(originCellIndex[czi][i]);
                            relayData[i + offset] = allTmpData[iter->second];
                        }
                        offset += originCellIndex[czi].size();
                    }
                }
                relayScatterFunc::relayScatter<real>(mesh().cellOfProc(), relayData, tmpData);
                relayScatterFunc::relayReorder<real>(tmpData, lastArrayArray_[i], mesh().perm(),
                                                     mesh().iperm());
            }
        }
    } else {
        realArray tmpData;
        realArray relayData;
        if (readLast || readFromGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.read(tmpData, groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_, "rhs0", groupName, this->name());
        } else {
            fos.read(tmpData, this->name());
            fos.readRealAttributeFromDataset(rhs0_, "rhs0", this->name());
        }
        const std::map<integer, integer> &indexMap = mesh().indexMap();
        const auto &originCellIndex = mesh().originCellIndex();
        relayData.resize(tmpData.size());

        integer offset = 0;

        for (integer czi = 0; czi < originCellIndex.size(); czi++) {
            for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                std::map<integer, integer>::const_iterator iter =
                    indexMap.find(originCellIndex[czi][i]);
                relayData[i + offset] = tmpData[iter->second];
            }
            offset += originCellIndex[czi].size();
        }

        relayScatterFunc::relayReorder<real>(relayData, *this, mesh().perm(), mesh().iperm());

        integer ng = 0;
        if (!fos.exist("nTimeGroups")) {
            fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
        }
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                realArray tmpData;
                string groupName = "timeGroupm";
                groupName += toString(i + 1);
                fos.read(tmpData, groupName, this->name());

                const std::map<integer, integer> &indexMap = mesh().indexMap();
                const auto &originCellIndex = mesh().originCellIndex();

                integer offset = 0;

                for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                    for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                        std::map<integer, integer>::const_iterator iter =
                            indexMap.find(originCellIndex[czi][i]);
                        relayData[i + offset] = tmpData[iter->second];
                    }
                    offset += originCellIndex[czi].size();
                }
                relayScatterFunc::relayReorder<real>(relayData, lastArrayArray_[i], mesh().perm(),
                                                     mesh().iperm());
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::real, OpenHurricane::cellMesh>::interpolateRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {
    if (HurMPI::parRun()) {
        realArray relayData;
        if (HurMPI::master()) {
            if (readLast || readFromGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.read(relayData, groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_, "rhs0", groupName, this->name());
            } else {
                fos.read(relayData, this->name());
                fos.readRealAttributeFromDataset(rhs0_, "rhs0", this->name());
            }
        }

        realArray tmpData;
        relayScatterFunc::relayScatter<real>(mesh().tarOfProc(), relayData, tmpData);

        const auto &nbr = mesh().sorKNN();
        const auto &sor = mesh().sorCellCentre();
        const auto &tar = mesh().cellCentre();

        meshInterpolation itp(nbr, sor, tar);
        itp.interpolate<real>(tmpData, *this);

        integer ng = 0;
        if (HurMPI::master()) {
            if (!fos.exist("nTimeGroups")) {
                fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
            }
        }
        HurMPI::bcast(&ng, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (readLast && ng > 0) {

            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(i + 1);
                    fos.read(relayData, groupName, this->name());
                }
                relayScatterFunc::relayScatter<real>(mesh().tarOfProc(), relayData, tmpData);
                itp.interpolate<real>(tmpData, lastArrayArray_[i]);
            }
        }
    } else {

        realArray tmpData;
        realArray relayData;
        if (readLast || readFromGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.read(relayData, groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_, "rhs0", groupName, this->name());
        } else {
            fos.read(relayData, this->name());
            fos.readRealAttributeFromDataset(rhs0_, "rhs0", this->name());
        }

        relayScatterFunc::relayScatter<real>(mesh().tarOfProc(), relayData, tmpData);

        const auto &nbr = mesh().sorKNN();
        const auto &sor = mesh().sorCellCentre();
        const auto &tar = mesh().cellCentre();

        meshInterpolation itp(nbr, sor, tar);
        itp.interpolate<real>(tmpData, *this);

        integer ng = 0;
        if (!fos.exist("nTimeGroups")) {
            fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
        }
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                string groupName = "timeGroupm";
                groupName += toString(i + 1);
                fos.read(relayData, groupName, this->name());

                relayScatterFunc::relayScatter<real>(mesh().tarOfProc(), relayData, tmpData);
                itp.interpolate<real>(tmpData, lastArrayArray_[i]);
            }
        }
    }
}


namespace OpenHurricane {
    template <> void geometryArray<vector2D, cellMesh>::writeOutput(fileOsstream &fos) const {
        if (argParse::checkWriteSolution()) {
            bool isok = true;

            for (integer n = 0; n < this->internalArraySize(); ++n) {
                if (isnan(this->operator[](n).x()) || isinf(this->operator[](n).x()) ||
                    isnan(this->operator[](n).y()) || isinf(this->operator[](n).y())) {
                    isok = false;
                    break;
                }
            }
            if (!isok) {
                LFatal("Not a number or infinity of %s", this->name().c_str());
            }
        }
        // Only write the internale filed's value.
        Array<vector2D>::writeToStream(fos, this->internalArraySize());
    }

    template <>
    void
    OpenHurricane::geometryArray<vector2D, cellMesh>::writeOutputByMaster(fileOsstream &fos) const {
        if (HurMPI::parRun()) {
            integerList nSizeL;
            integerList displs;
            integer allSize = 0;
            realArray rootF;
            realArray cmpF;

            nSizeL.resize(HurMPI::getProcSize(), Zero);
            nSizeL[HurMPI::getProcRank()] = this->internalArraySize();
            if (HurMPI::master()) {
                displs.resize(HurMPI::getProcSize(), Zero);
            }
            HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());

            if (HurMPI::master()) {
                for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                    displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
                }
                for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                    allSize += nSizeL[ip];
                }
            }
            cmpF.resize(this->internalArraySize());
            if (HurMPI::master()) {
                rootF.resize(allSize);
            }
            if (argParse::checkWriteSolution()) {
                bool isok = true;

                for (integer n = 0; n < this->internalArraySize(); ++n) {
                    if (isnan(this->operator[](n).x()) || isinf(this->operator[](n).x()) ||
                        isnan(this->operator[](n).y()) || isinf(this->operator[](n).y())) {
                        isok = false;
                        break;
                    }
                }
                if (!isok) {
                    LFatal("Not a number or infinity of %s", this->name().c_str());
                }
            }
            HurMPI::barrier(HurMPI::getComm());

            for (integer i = 0; i < feature<vector2D>::nElements_; ++i) {
                for (integer j = 0; j < this->internalArraySize(); ++j) {
                    cmpF[j] = this->operator[](j)[i];
                }

                HurMPI::Request request;
                HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                                 rootF.data(), nSizeL.data(), displs.data(), feature<real>::MPIType,
                                 HurMPI::masterNo(), HurMPI::getComm(), &request);
                HurMPI::wait(&request, MPI_STATUSES_IGNORE);

                if (HurMPI::master()) {
                    fos.write(reinterpret_cast<const char *>(&rootF[0]), rootF.byteSize());
                }
            }
        } else {
            writeOutput(fos);
        }
    }

    template <>
    void OpenHurricane::geometryArray<vector2D, cellMesh>::writeOutput(fileOsstream &fos,
                                                                   const integer fzid) const {
        const auto &displs = mesh().globalFaceZoneInfo(fzid).faceDispls();
        const auto &recvcnt = mesh().globalFaceZoneInfo(fzid).faceRecvcnt();
        const auto &faces = mesh().faces();
        const auto &fw = mesh().faceWgt();
        const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();
        const auto fsize = mesh().globalFaceZoneInfo(fzid).fZ().size();
        const auto nTotalFaces = mesh().globalFaceZoneInfo(fzid).totalFaces();

        realArray faceValue(fsize, Zero);
        realArray writeField;
        if (HurMPI::parRun()) {
            if (HurMPI::master()) {
                writeField.resize(nTotalFaces, Zero);
            }
        }
        for (int j = 0; j < vector2D::nElements_; ++j) {
            integer count = 0;
            for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
                const integer cl = faces[fi].leftCell();
                const integer cr = faces[fi].rightCell();
                faceValue[count++] =
                    fw[fi] * this->operator[](cl)[j] + (1.0 - fw[fi]) * this->operator[](cr)[j];
            }
            if (!HurMPI::parRun()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < fsize) {
                        sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                             << " ";
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (fsize) {
                        fos.write(reinterpret_cast<const char *>(&faceValue[0]),
                                  fsize * sizeof(real));
                    }
                }
                continue;
            }
            HurMPI::Request request;
            HurMPI::igatherv(faceValue.data(), fsize, feature<real>::MPIType, writeField.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeField.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                             << " ";
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeField.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeField[0]),
                                  writeField.size() * sizeof(real));
                    }
                }
            }
        }
    }

} // namespace OpenHurricane

template <>
void OpenHurricane::geometryArray<OpenHurricane::vector2D, OpenHurricane::cellMesh>::writeMinMaxOutput(
    fileOsstream &fos, const integer fzid) const {
    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();

    for (integer j = 0; j < vector2D::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
            const integer cl = faces[fi].leftCell();
            const integer cr = faces[fi].rightCell();
            real vv = 0.5 * this->operator[](cl)[j] + 0.5 * this->operator[](cr)[j];
            minX = min(minX, vv);
            maxX = max(maxX, vv);
        }
        if (HurMPI::parRun()) {
            HurMPI::reduce(minX, MPI_MIN);
            HurMPI::reduce(maxX, MPI_MAX);
        }
        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::vector2D, OpenHurricane::cellMesh>::writeRelay(
    hdf5O &fos, const bool writeLast, const bool writeToGroup) const {
    integerList nSizeL;
    integerList displs;
    integer allSize = 0;
    vector2DArray rootV2D;
    realArray rootF;
    realArray cmpF;
    if (HurMPI::parRun()) {
        nSizeL.resize(HurMPI::getProcSize(), Zero);
        nSizeL[HurMPI::getProcRank()] = this->internalArraySize();
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
        cmpF.resize(this->internalArraySize());
        if (HurMPI::master()) {
            rootF.resize(allSize);
            rootV2D.resize(allSize);
        }

        HurMPI::barrier(HurMPI::getComm());

        for (integer i = 0; i < feature<vector2D>::nElements_; ++i) {
            for (integer j = 0; j < this->internalArraySize(); ++j) {
                cmpF[j] = this->operator[](j)[i];
            }

            HurMPI::Request request;
            HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                             rootF.data(), nSizeL.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            for (integer j = 0; j < rootF.size(); ++j) {
                rootV2D[j][i] = rootF[j];
            }
        }
        if (HurMPI::master()) {
            if (writeLast || writeToGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.write(rootV2D, rootV2D.size(), groupName, this->name());
                fos.writeRealAttributeToDataset(rhs0_[0], "rhs00", groupName, this->name());
                fos.writeRealAttributeToDataset(rhs0_[1], "rhs01", groupName, this->name());
            } else {
                fos.write(rootV2D, rootV2D.size(), this->name());
                fos.writeRealAttributeToDataset(rhs0_[0], "rhs00", this->name());
                fos.writeRealAttributeToDataset(rhs0_[1], "rhs01", this->name());
            }
        }
    } else {
        if (writeLast || writeToGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.write(*this, this->internalArraySize(), groupName, this->name());
            fos.writeRealAttributeToDataset(rhs0_[0], "rhs00", groupName, this->name());
            fos.writeRealAttributeToDataset(rhs0_[1], "rhs01", groupName, this->name());
        } else {
            fos.write(*this, this->internalArraySize(), this->name());
            fos.writeRealAttributeToDataset(rhs0_[0], "rhs00", this->name());
            fos.writeRealAttributeToDataset(rhs0_[1], "rhs01", this->name());
        }
    }
    if (writeLast) {
        if (lastArrayArray_.size() == 0) {
            return;
        }
        if (HurMPI::parRun()) {
            for (integer k = 0; k < lastArrayArray_.size(); ++k) {
                HurMPI::barrier(HurMPI::getComm());
                for (integer i = 0; i < feature<vector2D>::nElements_; ++i) {
                    for (integer j = 0; j < this->internalArraySize(); ++j) {
                        cmpF[j] = lastArrayArray_[k][j][i];
                    }

                    HurMPI::Request request;
                    HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                                     rootF.data(), nSizeL.data(), displs.data(),
                                     feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm(),
                                     &request);
                    HurMPI::wait(&request, MPI_STATUSES_IGNORE);

                    for (integer j = 0; j < rootF.size(); ++j) {
                        rootV2D[j][i] = rootF[j];
                    }
                }
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(k + 1);
                    fos.write(rootV2D, rootV2D.size(), groupName, this->name());
                }
            }
        } else {
            for (integer k = 0; k < lastArrayArray_.size(); ++k) {
                string groupName = "timeGroupm";
                groupName += toString(k + 1);
                fos.write(lastArrayArray_[k], this->internalArraySize(), groupName, this->name());
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::vector2D, OpenHurricane::cellMesh>::readRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {
    if (HurMPI::parRun()) {
        vector2DArray allTmpData;
        vector2DArray relayData;
        if (HurMPI::master()) {
            if (readLast || readFromGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.read(allTmpData, groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", groupName, this->name());
            } else {
                fos.read(allTmpData, this->name());
                fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", this->name());
                fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", this->name());
            }
            const std::map<integer, integer> &indexMap = mesh().indexMap();
            const integerArrayArray &originCellIndex = mesh().originCellIndex();
            relayData.resize(allTmpData.size());

            integer offset = 0;

            for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                    std::map<integer, integer>::const_iterator iter =
                        indexMap.find(originCellIndex[czi][i]);
                    relayData[i + offset] = allTmpData[iter->second];
                }
                offset += originCellIndex[czi].size();
            }
        }
        vector2DArray tmpData;
        relayScatterFunc::relayScatter<vector2D>(mesh().cellOfProc(), relayData, tmpData);
        relayScatterFunc::relayReorder<vector2D>(tmpData, *this, mesh().perm(), mesh().iperm());
        integer ng = 0;
        if (HurMPI::master()) {
            if (!fos.exist("nTimeGroups")) {
                fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
            }
        }
        HurMPI::bcast(&ng, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(i + 1);
                    fos.read(allTmpData, groupName, this->name());
                    const std::map<integer, integer> &indexMap = mesh().indexMap();
                    const integerArrayArray &originCellIndex = mesh().originCellIndex();

                    integer offset = 0;

                    for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                        for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                            std::map<integer, integer>::const_iterator iter =
                                indexMap.find(originCellIndex[czi][i]);
                            relayData[i + offset] = allTmpData[iter->second];
                        }
                        offset += originCellIndex[czi].size();
                    }
                }
                relayScatterFunc::relayScatter<vector2D>(mesh().cellOfProc(), relayData, tmpData);
                relayScatterFunc::relayReorder<vector2D>(tmpData, lastArrayArray_[i], mesh().perm(),
                                                         mesh().iperm());
            }
        }
    } else {
        vector2DArray tmpData;
        vector2DArray relayData;
        if (readLast || readFromGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.read(tmpData, groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", groupName, this->name());
        } else {
            fos.read(tmpData, this->name());
            fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", this->name());
            fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", this->name());
        }
        const std::map<integer, integer> &indexMap = mesh().indexMap();
        const integerArrayArray &originCellIndex = mesh().originCellIndex();
        relayData.resize(tmpData.size());

        integer offset = 0;

        for (integer czi = 0; czi < originCellIndex.size(); czi++) {
            for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                std::map<integer, integer>::const_iterator iter =
                    indexMap.find(originCellIndex[czi][i]);
                relayData[i + offset] = tmpData[iter->second];
            }
            offset += originCellIndex[czi].size();
        }
        relayScatterFunc::relayReorder<vector2D>(relayData, *this, mesh().perm(), mesh().iperm());
        integer ng = 0;
        if (!fos.exist("nTimeGroups")) {
            fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
        }
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                vector2DArray tmpData;
                string groupName = "timeGroupm";
                groupName += toString(i + 1);
                fos.read(tmpData, groupName, this->name());
                const std::map<integer, integer> &indexMap = mesh().indexMap();
                const integerArrayArray &originCellIndex = mesh().originCellIndex();

                integer offset = 0;

                for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                    for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                        std::map<integer, integer>::const_iterator iter =
                            indexMap.find(originCellIndex[czi][i]);
                        relayData[i + offset] = tmpData[iter->second];
                    }
                    offset += originCellIndex[czi].size();
                }
                relayScatterFunc::relayReorder<vector2D>(relayData, lastArrayArray_[i],
                                                         mesh().perm(), mesh().iperm());
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::vector2D, OpenHurricane::cellMesh>::interpolateRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {
    if (HurMPI::parRun()) {
        vector2DArray relayData;
        if (HurMPI::master()) {
            if (readLast || readFromGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.read(relayData, groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", groupName, this->name());
            } else {
                fos.read(relayData, this->name());
                fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", this->name());
                fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", this->name());
            }
        }
        vector2DArray tmpData;
        relayScatterFunc::relayScatter<vector2D>(mesh().tarOfProc(), relayData, tmpData);

        const auto &nbr = mesh().sorKNN();
        const auto &sor = mesh().sorCellCentre();
        const auto &tar = mesh().cellCentre();

        meshInterpolation itp(nbr, sor, tar);
        itp.interpolate<vector2D>(tmpData, *this);

        integer ng = 0;
        if (HurMPI::master()) {
            if (!fos.exist("nTimeGroups")) {
                fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
            }
        }
        HurMPI::bcast(&ng, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(i + 1);
                    fos.read(relayData, groupName, this->name());
                }
                relayScatterFunc::relayScatter<vector2D>(mesh().tarOfProc(), relayData, tmpData);
                itp.interpolate<vector2D>(tmpData, lastArrayArray_[i]);
            }
        }
    } else {
        vector2DArray tmpData;
        vector2DArray relayData;
        if (readLast || readFromGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.read(relayData, groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", groupName, this->name());
        } else {
            fos.read(relayData, this->name());
            fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", this->name());
            fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", this->name());
        }
        relayScatterFunc::relayScatter<vector2D>(mesh().tarOfProc(), relayData, tmpData);

        const auto &nbr = mesh().sorKNN();
        const auto &sor = mesh().sorCellCentre();
        const auto &tar = mesh().cellCentre();

        meshInterpolation itp(nbr, sor, tar);
        itp.interpolate<vector2D>(tmpData, *this);

        integer ng = 0;
        if (!fos.exist("nTimeGroups")) {
            fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
        }
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                string groupName = "timeGroupm";
                groupName += toString(i + 1);
                fos.read(relayData, groupName, this->name());

                relayScatterFunc::relayScatter<vector2D>(mesh().tarOfProc(), relayData, tmpData);
                itp.interpolate<vector2D>(tmpData, lastArrayArray_[i]);
            }
        }
    }
}

namespace OpenHurricane {
    template <> void geometryArray<vector, cellMesh>::getAveAndMaxRHS() const {
        if (!rhsAvePtr_) {
            rhsAvePtr_.reset(new vector());
        }
        if (!rhsMaxPtr_) {
            rhsMaxPtr_.reset(new vector());
        }

        vector rhsAve = 0.0;
        vector rhsMax = 0.0;
        vector zRhs = 0.0;

        for (integer i = 0; i < mesh().nCells(); ++i) {
            zRhs[0] = fabs(rhs()[i][0]) / mesh().cellVol()[i];
            rhsAve[0] += sqr(zRhs[0]);
            rhsMax[0] = max(rhsMax[0], zRhs[0]);
            zRhs[1] = fabs(rhs()[i][1]) / mesh().cellVol()[i];
            rhsAve[1] += sqr(zRhs[1]);
            rhsMax[1] = max(rhsMax[1], zRhs[1]);
            zRhs[2] = fabs(rhs()[i][2]) / mesh().cellVol()[i];
            rhsAve[2] += sqr(zRhs[2]);
            rhsMax[2] = max(rhsMax[2], zRhs[2]);
        }

        realArray rhsAveL0(HurMPI::getProcSize(), Zero);
        HurMPI::gather(&rhsAve[0], 1, feature<real>::MPIType, rhsAveL0.data(), 1,
                       feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        realArray rhsAveL1(HurMPI::getProcSize(), Zero);
        HurMPI::gather(&rhsAve[1], 1, feature<real>::MPIType, rhsAveL1.data(), 1,
                       feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

        realArray rhsAveL2(HurMPI::getProcSize(), Zero);
        HurMPI::gather(&rhsAve[2], 1, feature<real>::MPIType, rhsAveL2.data(), 1,
                       feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm());

        if (HurMPI::master()) {
            rhsAve = 0.0;

            for (integer i = 0; i < rhsAveL0.size(); i++) {
                rhsAve[0] += rhsAveL0[i];
                rhsAve[1] += rhsAveL1[i];
                rhsAve[2] += rhsAveL2[i];
            }
            rhsAve /= real(mesh().allMeshCellNumber());
        }

        HurMPI::bcast(&rhsAve[0], 1, feature<real>::MPIType, HurMPI::masterNo());
        HurMPI::bcast(&rhsAve[1], 1, feature<real>::MPIType, HurMPI::masterNo());
        HurMPI::bcast(&rhsAve[2], 1, feature<real>::MPIType, HurMPI::masterNo());
        *rhsAvePtr_ = rhsAve;
        HurMPI::allReduceVS(rhsMax, MPI_MAX);
        *rhsMaxPtr_ = rhsMax;
    }

    template <> void geometryArray<vector, cellMesh>::updateBoundary(integer layerI) {
        vectorTransfer myTransfer(this->mesh(), *this, false, false);
        myTransfer.transferInit(layerI);
        if (!boundariesPtr_) {
            const controller &cont = tb().Iteration().cont();
            boundariesPtr_ =
                new boundaryCondition<vector, cellMesh>(mesh_.faceZones(), *this, cont);
        }
        boundariesPtr_->calculate();
        myTransfer.transferring(layerI);
    }

    template <> vector geometryArray<vector, cellMesh>::initialize() {
        if (!boundariesPtr_) {
            const controller &cont = tb().Iteration().cont();
            boundariesPtr_ =
                new boundaryCondition<vector, cellMesh>(mesh_.faceZones(), *this, cont);
        }
        const controller &cont = tb().Iteration().cont();
        if (!cont.found("initialization")) {
            LFatal("Cannot find initialization controller in: %s", cont.name().c_str());
        }
        const auto &initCont = cont.subController("initialization");
        if (initCont.found("initFromValue")) {
            const auto &initValueCont = initCont.subController("initFromValue");

            if (initValueCont.found(this->name())) {
                const vector initValue = initValueCont.findType<vector>(this->name(), vector(0.0));
                *this = initValue;
                return initValue;
            }
        }
        for (integer i = 0; i < boundariesPtr_->size(); i++) {
            if ((*boundariesPtr_)[i].boundaryfZ().name() ==
                cont.subController("initialization").findWord("initFromBoundary")) {
                *this = ((*boundariesPtr_)[i].initValue());
                return (*boundariesPtr_)[i].initValue();
            }
        }
        checkWarningStr(("The variable: " + this->name() + " is not initialized."));
        return vector(Zero);
    }

    template <>
    void OpenHurricane::geometryArray<vector, cellMesh>::writeOutput(fileOsstream &fos) const {
        if (argParse::checkWriteSolution()) {
            bool isok = true;

            for (integer n = 0; n < this->internalArraySize(); ++n) {
                if (isnan(this->operator[](n).x()) || isinf(this->operator[](n).x()) ||
                    isnan(this->operator[](n).y()) || isinf(this->operator[](n).y()) ||
                    isnan(this->operator[](n).z()) || isinf(this->operator[](n).z())) {
                    isok = false;
                    break;
                }
            }
            if (!isok) {
                LFatal("Not a number or infinity of %s", this->name().c_str());
            }
        }
        // Only write the internale filed's value.
        Array<vector>::writeToStream(fos, this->internalArraySize());
    }

    template <>
    void OpenHurricane::geometryArray<vector, cellMesh>::writeOutputByMaster(fileOsstream &fos) const {
        if (HurMPI::parRun()) {
            integerList nSizeL;
            integerList displs;
            integer allSize = 0;
            realArray rootF;
            realArray cmpF;

            nSizeL.resize(HurMPI::getProcSize(), Zero);
            nSizeL[HurMPI::getProcRank()] = this->internalArraySize();
            if (HurMPI::master()) {
                displs.resize(HurMPI::getProcSize(), Zero);
            }
            HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());

            if (HurMPI::master()) {
                for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                    displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
                }
                for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                    allSize += nSizeL[ip];
                }
            }
            cmpF.resize(this->internalArraySize());
            if (HurMPI::master()) {
                rootF.resize(allSize);
            }
            if (argParse::checkWriteSolution()) {
                bool isok = true;

                for (integer n = 0; n < this->internalArraySize(); ++n) {
                    if (isnan(this->operator[](n).x()) || isinf(this->operator[](n).x()) ||
                        isnan(this->operator[](n).y()) || isinf(this->operator[](n).y()) ||
                        isnan(this->operator[](n).z()) || isinf(this->operator[](n).z())) {
                        isok = false;
                        break;
                    }
                }
                if (!isok) {
                    LFatal("Not a number or infinity of %s", this->name().c_str());
                }
            }

            HurMPI::barrier(HurMPI::getComm());

            for (integer i = 0; i < feature<vector>::nElements_; ++i) {
                for (integer j = 0; j < this->internalArraySize(); ++j) {
                    cmpF[j] = this->operator[](j)[i];
                }

                HurMPI::Request request;
                HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                                 rootF.data(), nSizeL.data(), displs.data(), feature<real>::MPIType,
                                 HurMPI::masterNo(), HurMPI::getComm(), &request);
                HurMPI::wait(&request, MPI_STATUSES_IGNORE);

                if (HurMPI::master()) {
                    fos.write(reinterpret_cast<const char *>(&rootF[0]), rootF.byteSize());
                }
            }
        } else {
            writeOutput(fos);
        }
    }

    template <>
    void OpenHurricane::geometryArray<vector, cellMesh>::writeOutput(fileOsstream &fos,
                                                                 const integer fzid) const {
        const auto &displs = mesh().globalFaceZoneInfo(fzid).faceDispls();
        const auto &recvcnt = mesh().globalFaceZoneInfo(fzid).faceRecvcnt();
        const auto &faces = mesh().faces();
        const auto &fw = mesh().faceWgt();
        const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();
        const auto fsize = mesh().globalFaceZoneInfo(fzid).fZ().size();
        const auto nTotalFaces = mesh().globalFaceZoneInfo(fzid).totalFaces();

        realArray faceValue(fsize, Zero);
        realArray writeField;
        if (HurMPI::parRun()) {
            if (HurMPI::master()) {
                writeField.resize(nTotalFaces, Zero);
            }
        }
        for (int j = 0; j < vector::nElements_; ++j) {
            integer count = 0;
            for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
                const integer cl = faces[fi].leftCell();
                const integer cr = faces[fi].rightCell();
                faceValue[count++] = 0.5 * this->operator[](cl)[j] + 0.5 * this->operator[](cr)[j];
            }
            if (!HurMPI::parRun()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < fsize) {
                        sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                             << " ";
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (fsize) {
                        fos.write(reinterpret_cast<const char *>(&faceValue[0]),
                                  fsize * sizeof(real));
                    }
                }
                continue;
            }

            HurMPI::Request request;
            HurMPI::igatherv(faceValue.data(), fsize, feature<real>::MPIType, writeField.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeField.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                             << " ";
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeField.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeField[0]),
                                  writeField.size() * sizeof(real));
                    }
                }
            }
        }
    }
} // namespace OpenHurricane

template <>
void OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::cellMesh>::writeMinMaxOutput(
    fileOsstream &fos, const integer fzid) const {
    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();

    for (integer j = 0; j < vector::nElements_; ++j) {

        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
            const integer cl = faces[fi].leftCell();
            const integer cr = faces[fi].rightCell();
            real vv = 0.5 * this->operator[](cl)[j] + 0.5 * this->operator[](cr)[j];
            minX = min(minX, vv);
            maxX = max(maxX, vv);
        }
        if (HurMPI::parRun()) {
            HurMPI::reduce(minX, MPI_MIN);
            HurMPI::reduce(maxX, MPI_MAX);
        }
        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::cellMesh>::writeRelay(
    hdf5O &fos, const bool writeLast, const bool writeToGroup) const {
    integerList nSizeL;
    integerList displs;
    integer allSize;
    vectorArray rootV;
    realArray rootF;
    realArray cmpF;
    if (HurMPI::parRun()) {
        nSizeL.resize(HurMPI::getProcSize(), Zero);
        nSizeL[HurMPI::getProcRank()] = this->internalArraySize();
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
        cmpF.resize(this->internalArraySize());
        if (HurMPI::master()) {
            rootF.resize(allSize);
            rootV.resize(allSize);
        }

        HurMPI::barrier(HurMPI::getComm());

        for (integer i = 0; i < feature<vector>::nElements_; ++i) {
            for (integer j = 0; j < this->internalArraySize(); ++j) {
                cmpF[j] = this->operator[](j)[i];
            }

            HurMPI::Request request;
            HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                             rootF.data(), nSizeL.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            for (integer j = 0; j < rootF.size(); ++j) {
                rootV[j][i] = rootF[j];
            }
        }
        if (HurMPI::master()) {
            if (writeLast || writeToGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.write(rootV, rootV.size(), groupName, this->name());
                fos.writeRealAttributeToDataset(rhs0_[0], "rhs00", groupName, this->name());
                fos.writeRealAttributeToDataset(rhs0_[1], "rhs01", groupName, this->name());
                fos.writeRealAttributeToDataset(rhs0_[2], "rhs02", groupName, this->name());
            } else {
                fos.write(rootV, rootV.size(), this->name());
                fos.writeRealAttributeToDataset(rhs0_[0], "rhs00", this->name());
                fos.writeRealAttributeToDataset(rhs0_[1], "rhs01", this->name());
                fos.writeRealAttributeToDataset(rhs0_[2], "rhs02", this->name());
            }
        }
    } else {
        if (writeLast || writeToGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.write(*this, this->internalArraySize(), groupName, this->name());
            fos.writeRealAttributeToDataset(rhs0_[0], "rhs00", groupName, this->name());
            fos.writeRealAttributeToDataset(rhs0_[1], "rhs01", groupName, this->name());
            fos.writeRealAttributeToDataset(rhs0_[2], "rhs02", groupName, this->name());
        } else {
            fos.write(*this, this->internalArraySize(), this->name());
            fos.writeRealAttributeToDataset(rhs0_[0], "rhs00", this->name());
            fos.writeRealAttributeToDataset(rhs0_[1], "rhs01", this->name());
            fos.writeRealAttributeToDataset(rhs0_[2], "rhs02", this->name());
        }
    }
    if (writeLast) {
        if (lastArrayArray_.size() == 0) {
            return;
        }
        if (HurMPI::parRun()) {
            for (integer k = 0; k < lastArrayArray_.size(); ++k) {
                HurMPI::barrier(HurMPI::getComm());
                for (integer i = 0; i < feature<vector>::nElements_; ++i) {
                    for (integer j = 0; j < this->internalArraySize(); ++j) {
                        cmpF[j] = lastArrayArray_[k][j][i];
                    }

                    HurMPI::Request request;
                    HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                                     rootF.data(), nSizeL.data(), displs.data(),
                                     feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm(),
                                     &request);
                    HurMPI::wait(&request, MPI_STATUSES_IGNORE);

                    for (integer j = 0; j < rootF.size(); ++j) {
                        rootV[j][i] = rootF[j];
                    }
                }
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(k + 1);
                    fos.write(rootV, rootV.size(), groupName, this->name());
                }
            }
        } else {
            for (integer k = 0; k < lastArrayArray_.size(); ++k) {
                string groupName = "timeGroupm";
                groupName += toString(k + 1);
                fos.write(lastArrayArray_[k], this->internalArraySize(), groupName, this->name());
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::cellMesh>::readRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {
    if (HurMPI::parRun()) {
        vectorArray allTmpData;
        vectorArray relayData;
        if (HurMPI::master()) {
            if (readLast || readFromGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.read(allTmpData, groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_[2], "rhs02", groupName, this->name());
            } else {
                fos.read(allTmpData, this->name());
                fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", this->name());
                fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", this->name());
                fos.readRealAttributeFromDataset(rhs0_[2], "rhs02", this->name());
            }
            const std::map<integer, integer> &indexMap = mesh().indexMap();
            const integerArrayArray &originCellIndex = mesh().originCellIndex();
            relayData.resize(allTmpData.size());

            integer offset = 0;

            for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                    std::map<integer, integer>::const_iterator iter =
                        indexMap.find(originCellIndex[czi][i]);
                    relayData[i + offset] = allTmpData[iter->second];
                }
                offset += originCellIndex[czi].size();
            }
        }
        vectorArray tmpData;
        relayScatterFunc::relayScatter<vector>(mesh().cellOfProc(), relayData, tmpData);
        relayScatterFunc::relayReorder<vector>(tmpData, *this, mesh().perm(), mesh().iperm());
        integer ng = 0;
        if (HurMPI::master()) {
            if (!fos.exist("nTimeGroups")) {
                fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
            }
        }
        HurMPI::bcast(&ng, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(i + 1);
                    fos.read(allTmpData, groupName, this->name());
                    const std::map<integer, integer> &indexMap = mesh().indexMap();
                    const integerArrayArray &originCellIndex = mesh().originCellIndex();

                    integer offset = 0;

                    for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                        for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                            std::map<integer, integer>::const_iterator iter =
                                indexMap.find(originCellIndex[czi][i]);
                            relayData[i + offset] = allTmpData[iter->second];
                        }
                        offset += originCellIndex[czi].size();
                    }
                }
                relayScatterFunc::relayScatter<vector>(mesh().cellOfProc(), relayData, tmpData);
                relayScatterFunc::relayReorder<vector>(tmpData, lastArrayArray_[i], mesh().perm(),
                                                       mesh().iperm());
            }
        }
    } else {
        vectorArray tmpData;
        vectorArray relayData;
        if (readLast || readFromGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.read(tmpData, groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_[2], "rhs02", groupName, this->name());
        } else {
            fos.read(tmpData, this->name());
            fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", this->name());
            fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", this->name());
            fos.readRealAttributeFromDataset(rhs0_[2], "rhs02", this->name());
        }
        const std::map<integer, integer> &indexMap = mesh().indexMap();
        const integerArrayArray &originCellIndex = mesh().originCellIndex();
        relayData.resize(tmpData.size());

        integer offset = 0;

        for (integer czi = 0; czi < originCellIndex.size(); czi++) {
            for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                std::map<integer, integer>::const_iterator iter =
                    indexMap.find(originCellIndex[czi][i]);
                relayData[i + offset] = tmpData[iter->second];
            }
            offset += originCellIndex[czi].size();
        }
        relayScatterFunc::relayReorder<vector>(relayData, *this, mesh().perm(), mesh().iperm());
        integer ng = 0;
        if (!fos.exist("nTimeGroups")) {
            fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
        }
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                vectorArray tmpData;
                string groupName = "timeGroupm";
                groupName += toString(i + 1);
                fos.read(tmpData, groupName, this->name());
                const std::map<integer, integer> &indexMap = mesh().indexMap();
                const integerArrayArray &originCellIndex = mesh().originCellIndex();

                integer offset = 0;

                for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                    for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                        std::map<integer, integer>::const_iterator iter =
                            indexMap.find(originCellIndex[czi][i]);
                        relayData[i + offset] = tmpData[iter->second];
                    }
                    offset += originCellIndex[czi].size();
                }
                relayScatterFunc::relayReorder<vector>(relayData, lastArrayArray_[i], mesh().perm(),
                                                       mesh().iperm());
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::vector, OpenHurricane::cellMesh>::interpolateRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {
    if (HurMPI::parRun()) {
        vectorArray relayData;
        if (HurMPI::master()) {
            if (readLast || readFromGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.read(relayData, groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", groupName, this->name());
                fos.readRealAttributeFromDataset(rhs0_[2], "rhs02", groupName, this->name());
            } else {
                fos.read(relayData, this->name());
                fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", this->name());
                fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", this->name());
                fos.readRealAttributeFromDataset(rhs0_[2], "rhs02", this->name());
            }
        }
        vectorArray tmpData;
        relayScatterFunc::relayScatter<vector>(mesh().tarOfProc(), relayData, tmpData);

        const auto &nbr = mesh().sorKNN();
        const auto &sor = mesh().sorCellCentre();
        const auto &tar = mesh().cellCentre();

        meshInterpolation itp(nbr, sor, tar);
        itp.interpolate<vector>(tmpData, *this);

        integer ng = 0;
        if (HurMPI::master()) {
            if (!fos.exist("nTimeGroups")) {
                fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
            }
        }
        HurMPI::bcast(&ng, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(i + 1);
                    fos.read(relayData, groupName, this->name());
                }

                relayScatterFunc::relayScatter<vector>(mesh().tarOfProc(), relayData, tmpData);
                itp.interpolate<vector>(tmpData, lastArrayArray_[i]);
            }
        }
    } else {
        vectorArray tmpData;
        vectorArray relayData;
        if (readLast || readFromGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.read(relayData, groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", groupName, this->name());
            fos.readRealAttributeFromDataset(rhs0_[2], "rhs02", groupName, this->name());
        } else {
            fos.read(relayData, this->name());
            fos.readRealAttributeFromDataset(rhs0_[0], "rhs00", this->name());
            fos.readRealAttributeFromDataset(rhs0_[1], "rhs01", this->name());
            fos.readRealAttributeFromDataset(rhs0_[2], "rhs02", this->name());
        }

        relayScatterFunc::relayScatter<vector>(mesh().tarOfProc(), relayData, tmpData);

        const auto &nbr = mesh().sorKNN();
        const auto &sor = mesh().sorCellCentre();
        const auto &tar = mesh().cellCentre();

        meshInterpolation itp(nbr, sor, tar);
        itp.interpolate<vector>(tmpData, *this);

        integer ng = 0;
        if (!fos.exist("nTimeGroups")) {
            fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
        }
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                string groupName = "timeGroupm";
                groupName += toString(i + 1);
                fos.read(relayData, groupName, this->name());

                relayScatterFunc::relayScatter<vector>(mesh().tarOfProc(), relayData, tmpData);
                itp.interpolate<vector>(tmpData, lastArrayArray_[i]);
            }
        }
    }
}

namespace OpenHurricane {
    template <> void geometryArray<tensor, cellMesh>::writeOutput(fileOsstream &fos) const {
        if (argParse::checkWriteSolution()) {
            bool isok = true;

            for (integer n = 0; n < this->internalArraySize(); ++n) {
                if (isnan(this->operator[](n).xx()) || isinf(this->operator[](n).xx()) ||
                    isnan(this->operator[](n).xy()) || isinf(this->operator[](n).xy()) ||
                    isnan(this->operator[](n).xz()) || isinf(this->operator[](n).xz()) ||
                    isnan(this->operator[](n).yx()) || isinf(this->operator[](n).yx()) ||
                    isnan(this->operator[](n).yy()) || isinf(this->operator[](n).yy()) ||
                    isnan(this->operator[](n).yz()) || isinf(this->operator[](n).yz()) ||
                    isnan(this->operator[](n).zx()) || isinf(this->operator[](n).zx()) ||
                    isnan(this->operator[](n).zy()) || isinf(this->operator[](n).zy()) ||
                    isnan(this->operator[](n).zz()) || isinf(this->operator[](n).zz())) {
                    isok = false;
                    break;
                }
            }
            if (!isok) {
                LFatal("Not a number or infinity of %s", this->name().c_str());
            }
        }
        // Only write the internale filed's value.
        Array<tensor>::writeToStream(fos, this->internalArraySize());
    }

    template <>
    void OpenHurricane::geometryArray<tensor, cellMesh>::writeOutputByMaster(fileOsstream &fos) const {
        if (HurMPI::parRun()) {
            integerList nSizeL;
            integerList displs;
            integer allSize = 0;
            realArray rootF;
            realArray cmpF;

            nSizeL.resize(HurMPI::getProcSize(), Zero);
            nSizeL[HurMPI::getProcRank()] = this->internalArraySize();
            if (HurMPI::master()) {
                displs.resize(HurMPI::getProcSize(), Zero);
            }
            HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());

            if (HurMPI::master()) {
                for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                    displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
                }
                for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                    allSize += nSizeL[ip];
                }
            }
            cmpF.resize(this->internalArraySize());
            if (HurMPI::master()) {
                rootF.resize(allSize);
            }
            if (argParse::checkWriteSolution()) {
                bool isok = true;

                for (integer n = 0; n < this->internalArraySize(); ++n) {
                    if (isnan(this->operator[](n).xx()) || isinf(this->operator[](n).xx()) ||
                        isnan(this->operator[](n).xy()) || isinf(this->operator[](n).xy()) ||
                        isnan(this->operator[](n).xz()) || isinf(this->operator[](n).xz()) ||
                        isnan(this->operator[](n).yx()) || isinf(this->operator[](n).yx()) ||
                        isnan(this->operator[](n).yy()) || isinf(this->operator[](n).yy()) ||
                        isnan(this->operator[](n).yz()) || isinf(this->operator[](n).yz()) ||
                        isnan(this->operator[](n).zx()) || isinf(this->operator[](n).zx()) ||
                        isnan(this->operator[](n).zy()) || isinf(this->operator[](n).zy()) ||
                        isnan(this->operator[](n).zz()) || isinf(this->operator[](n).zz())) {
                        isok = false;
                        break;
                    }
                }
                if (!isok) {
                    LFatal("Not a number or infinity of %s", this->name().c_str());
                }
            }

            HurMPI::barrier(HurMPI::getComm());

            for (integer i = 0; i < feature<tensor>::nElements_; ++i) {
                for (integer j = 0; j < this->internalArraySize(); ++j) {
                    cmpF[j] = this->operator[](j)[i];
                }

                HurMPI::Request request;
                HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                                 rootF.data(), nSizeL.data(), displs.data(), feature<real>::MPIType,
                                 HurMPI::masterNo(), HurMPI::getComm(), &request);
                HurMPI::wait(&request, MPI_STATUSES_IGNORE);

                if (HurMPI::master()) {
                    fos.write(reinterpret_cast<const char *>(&rootF[0]), rootF.byteSize());
                }
            }
        } else {
            writeOutput(fos);
        }
    }

    template <>
    void OpenHurricane::geometryArray<tensor, cellMesh>::writeOutput(fileOsstream &fos,
                                                                 const integer fzid) const {
        const auto &displs = mesh().globalFaceZoneInfo(fzid).faceDispls();
        const auto &recvcnt = mesh().globalFaceZoneInfo(fzid).faceRecvcnt();
        const auto &faces = mesh().faces();
        const auto &fw = mesh().faceWgt();
        const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();
        const auto fsize = mesh().globalFaceZoneInfo(fzid).fZ().size();
        const auto nTotalFaces = mesh().globalFaceZoneInfo(fzid).totalFaces();

        realArray faceValue(fsize, Zero);
        realArray writeField;
        if (HurMPI::parRun()) {
            if (HurMPI::master()) {
                writeField.resize(nTotalFaces, Zero);
            }
        }
        for (int j = 0; j < tensor::nElements_; ++j) {
            integer count = 0;
            for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
                const integer cl = faces[fi].leftCell();
                const integer cr = faces[fi].rightCell();
                faceValue[count++] =
                    fw[fi] * this->operator[](cl)[j] + (1.0 - fw[fi]) * this->operator[](cr)[j];
            }
            if (!HurMPI::parRun()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < fsize) {
                        sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                             << " ";
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (fsize) {
                        fos.write(reinterpret_cast<const char *>(&faceValue[0]),
                                  fsize * sizeof(real));
                    }
                }
                continue;
            }

            HurMPI::Request request;
            HurMPI::igatherv(faceValue.data(), fsize, feature<real>::MPIType, writeField.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeField.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                             << " ";
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeField.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeField[0]),
                                  writeField.size() * sizeof(real));
                    }
                }
            }
        }
    }

} // namespace OpenHurricane

template <>
void OpenHurricane::geometryArray<OpenHurricane::tensor, OpenHurricane::cellMesh>::writeMinMaxOutput(
    fileOsstream &fos, const integer fzid) const {
    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();

    for (integer j = 0; j < tensor::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
            const integer cl = faces[fi].leftCell();
            const integer cr = faces[fi].rightCell();
            real vv = 0.5 * this->operator[](cl)[j] + 0.5 * this->operator[](cr)[j];

            minX = min(minX, vv);
            maxX = max(maxX, vv);
        }
        if (HurMPI::parRun()) {
            HurMPI::reduce(minX, MPI_MIN);
            HurMPI::reduce(maxX, MPI_MAX);
        }
        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::tensor, OpenHurricane::cellMesh>::writeRelay(
    hdf5O &fos, const bool writeLast, const bool writeToGroup) const {
    integerList nSizeL;
    integerList displs;
    integer allSize;
    tensorArray rootV;
    realArray rootF;
    realArray cmpF;
    if (HurMPI::parRun()) {
        nSizeL.resize(HurMPI::getProcSize(), Zero);
        nSizeL[HurMPI::getProcRank()] = this->internalArraySize();
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
        cmpF.resize(this->internalArraySize());
        if (HurMPI::master()) {
            rootF.resize(allSize);
            rootV.resize(allSize);
        }

        HurMPI::barrier(HurMPI::getComm());

        for (integer i = 0; i < feature<tensor>::nElements_; ++i) {
            for (integer j = 0; j < this->internalArraySize(); ++j) {
                cmpF[j] = this->operator[](j)[i];
            }

            HurMPI::Request request;
            HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                             rootF.data(), nSizeL.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            for (integer j = 0; j < rootF.size(); ++j) {
                rootV[j][i] = rootF[j];
            }
        }
        if (HurMPI::master()) {
            if (writeLast || writeToGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.write(rootV, real(1), rootV.size(), groupName, this->name());
                for (integer i = 0; i < feature<tensor>::nElements_; ++i) {
                    string rhsName = "rhs0" + toString(i);
                    fos.writeRealAttributeToDataset(rhs0_[i], rhsName, groupName, this->name());
                }
            } else {
                fos.write(rootV, real(1), rootV.size(), this->name());
                for (integer i = 0; i < feature<tensor>::nElements_; ++i) {
                    string rhsName = "rhs0" + toString(i);
                    fos.writeRealAttributeToDataset(rhs0_[i], rhsName, this->name());
                }
            }
        }
    } else {
        if (writeLast || writeToGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.write(*this, real(1), this->internalArraySize(), groupName, this->name());
            for (integer i = 0; i < feature<tensor>::nElements_; ++i) {
                string rhsName = "rhs0" + toString(i);
                fos.writeRealAttributeToDataset(rhs0_[i], rhsName, groupName, this->name());
            }
        } else {
            fos.write(*this, real(1), this->internalArraySize(), this->name());
            for (integer i = 0; i < feature<tensor>::nElements_; ++i) {
                string rhsName = "rhs0" + toString(i);
                fos.writeRealAttributeToDataset(rhs0_[i], rhsName, this->name());
            }
        }
    }
    if (writeLast) {
        if (lastArrayArray_.size() == 0) {
            return;
        }
        if (HurMPI::parRun()) {
            for (integer k = 0; k < lastArrayArray_.size(); ++k) {
                HurMPI::barrier(HurMPI::getComm());
                for (integer i = 0; i < feature<tensor>::nElements_; ++i) {
                    for (integer j = 0; j < this->internalArraySize(); ++j) {
                        cmpF[j] = lastArrayArray_[k][j][i];
                    }

                    HurMPI::Request request;
                    HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                                     rootF.data(), nSizeL.data(), displs.data(),
                                     feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm(),
                                     &request);
                    HurMPI::wait(&request, MPI_STATUSES_IGNORE);

                    for (integer j = 0; j < rootF.size(); ++j) {
                        rootV[j][i] = rootF[j];
                    }
                }
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(k + 1);
                    fos.write(rootV, rootV.size(), groupName, this->name());
                }
            }
        } else {
            for (integer k = 0; k < lastArrayArray_.size(); ++k) {
                string groupName = "timeGroupm";
                groupName += toString(k + 1);
                fos.write(lastArrayArray_[k], this->internalArraySize(), groupName, this->name());
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::tensor, OpenHurricane::cellMesh>::readRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {
    if (HurMPI::parRun()) {
        tensorArray allTmpData;
        tensorArray relayData;
        if (HurMPI::master()) {
            if (readLast || readFromGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.read(allTmpData, groupName, this->name());
                for (integer i = 0; i < feature<tensor>::nElements_; i++) {
                    string rhsName = "rhs0" + toString(i);
                    fos.readRealAttributeFromDataset(rhs0_[i], rhsName, groupName, this->name());
                }
            } else {
                fos.read(allTmpData, this->name());
                for (integer i = 0; i < feature<tensor>::nElements_; i++) {
                    string rhsName = "rhs0" + toString(i);
                    fos.readRealAttributeFromDataset(rhs0_[i], rhsName, this->name());
                }
            }
            const std::map<integer, integer> &indexMap = mesh().indexMap();
            const integerArrayArray &originCellIndex = mesh().originCellIndex();
            relayData.resize(allTmpData.size());

            integer offset = 0;

            for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                    std::map<integer, integer>::const_iterator iter =
                        indexMap.find(originCellIndex[czi][i]);
                    relayData[i + offset] = allTmpData[iter->second];
                }
                offset += originCellIndex[czi].size();
            }
        }
        tensorArray tmpData;
        relayScatterFunc::relayScatter<tensor>(mesh().cellOfProc(), relayData, tmpData);
        relayScatterFunc::relayReorder<tensor>(tmpData, *this, mesh().perm(), mesh().iperm());
        integer ng = 0;
        if (HurMPI::master()) {
            if (!fos.exist("nTimeGroups")) {
                fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
            }
        }
        HurMPI::bcast(&ng, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(i + 1);
                    fos.read(allTmpData, groupName, this->name());
                    const std::map<integer, integer> &indexMap = mesh().indexMap();
                    const integerArrayArray &originCellIndex = mesh().originCellIndex();

                    integer offset = 0;

                    for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                        for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                            std::map<integer, integer>::const_iterator iter =
                                indexMap.find(originCellIndex[czi][i]);
                            relayData[i + offset] = allTmpData[iter->second];
                        }
                        offset += originCellIndex[czi].size();
                    }
                }
                relayScatterFunc::relayScatter<tensor>(mesh().cellOfProc(), relayData, tmpData);
                relayScatterFunc::relayReorder<tensor>(tmpData, lastArrayArray_[i], mesh().perm(),
                                                       mesh().iperm());
            }
        }
    } else {
        tensorArray tmpData;
        tensorArray relayData;
        if (readLast || readFromGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.read(tmpData, groupName, this->name());
            for (integer i = 0; i < feature<tensor>::nElements_; i++) {
                string rhsName = "rhs0" + toString(i);
                fos.readRealAttributeFromDataset(rhs0_[i], rhsName, groupName, this->name());
            }
        } else {
            fos.read(tmpData, this->name());
            for (integer i = 0; i < feature<tensor>::nElements_; i++) {
                string rhsName = "rhs0" + toString(i);
                fos.readRealAttributeFromDataset(rhs0_[i], rhsName, this->name());
            }
        }
        const std::map<integer, integer> &indexMap = mesh().indexMap();
        const integerArrayArray &originCellIndex = mesh().originCellIndex();
        relayData.resize(tmpData.size());

        integer offset = 0;

        for (integer czi = 0; czi < originCellIndex.size(); czi++) {
            for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                std::map<integer, integer>::const_iterator iter =
                    indexMap.find(originCellIndex[czi][i]);
                relayData[i + offset] = tmpData[iter->second];
            }
            offset += originCellIndex[czi].size();
        }
        relayScatterFunc::relayReorder<tensor>(relayData, *this, mesh().perm(), mesh().iperm());
        integer ng = 0;
        if (!fos.exist("nTimeGroups")) {
            fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
        }
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                tensorArray tmpData;
                string groupName = "timeGroupm";
                groupName += toString(i + 1);
                fos.read(tmpData, groupName, this->name());
                const std::map<integer, integer> &indexMap = mesh().indexMap();
                const integerArrayArray &originCellIndex = mesh().originCellIndex();

                integer offset = 0;

                for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                    for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                        std::map<integer, integer>::const_iterator iter =
                            indexMap.find(originCellIndex[czi][i]);
                        relayData[i + offset] = tmpData[iter->second];
                    }
                    offset += originCellIndex[czi].size();
                }
                relayScatterFunc::relayReorder<tensor>(relayData, lastArrayArray_[i], mesh().perm(),
                                                       mesh().iperm());
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::tensor, OpenHurricane::cellMesh>::interpolateRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {
    if (HurMPI::parRun()) {
        tensorArray relayData;
        if (HurMPI::master()) {
            if (readLast || readFromGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.read(relayData, groupName, this->name());
                for (integer i = 0; i < feature<tensor>::nElements_; i++) {
                    string rhsName = "rhs0" + toString(i);
                    fos.readRealAttributeFromDataset(rhs0_[i], rhsName, groupName, this->name());
                }
            } else {
                fos.read(relayData, this->name());
                for (integer i = 0; i < feature<tensor>::nElements_; i++) {
                    string rhsName = "rhs0" + toString(i);
                    fos.readRealAttributeFromDataset(rhs0_[i], rhsName, this->name());
                }
            }
        }
        tensorArray tmpData;
        relayScatterFunc::relayScatter<tensor>(mesh().tarOfProc(), relayData, tmpData);

        const auto &nbr = mesh().sorKNN();
        const auto &sor = mesh().sorCellCentre();
        const auto &tar = mesh().cellCentre();

        meshInterpolation itp(nbr, sor, tar);
        itp.interpolate<tensor>(tmpData, *this);

        integer ng = 0;
        if (HurMPI::master()) {
            if (!fos.exist("nTimeGroups")) {
                fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
            }
        }
        HurMPI::bcast(&ng, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(i + 1);
                    fos.read(relayData, groupName, this->name());
                }
                relayScatterFunc::relayScatter<tensor>(mesh().tarOfProc(), relayData, tmpData);
                itp.interpolate<tensor>(tmpData, lastArrayArray_[i]);
            }
        }
    } else {
        tensorArray tmpData;
        tensorArray relayData;
        if (readLast || readFromGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.read(relayData, groupName, this->name());
            for (integer i = 0; i < feature<tensor>::nElements_; i++) {
                string rhsName = "rhs0" + toString(i);
                fos.readRealAttributeFromDataset(rhs0_[i], rhsName, groupName, this->name());
            }
        } else {
            fos.read(relayData, this->name());
            for (integer i = 0; i < feature<tensor>::nElements_; i++) {
                string rhsName = "rhs0" + toString(i);
                fos.readRealAttributeFromDataset(rhs0_[i], rhsName, this->name());
            }
        }
        relayScatterFunc::relayReorder<tensor>(relayData, *this, mesh().perm(), mesh().iperm());

        const auto &nbr = mesh().sorKNN();
        const auto &sor = mesh().sorCellCentre();
        const auto &tar = mesh().cellCentre();

        meshInterpolation itp(nbr, sor, tar);
        itp.interpolate<tensor>(tmpData, *this);

        integer ng = 0;
        if (!fos.exist("nTimeGroups")) {
            fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
        }
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                string groupName = "timeGroupm";
                groupName += toString(i + 1);
                fos.read(relayData, groupName, this->name());

                relayScatterFunc::relayReorder<tensor>(relayData, *this, mesh().perm(),
                                                       mesh().iperm());
                itp.interpolate<tensor>(tmpData, lastArrayArray_[i]);
            }
        }
    }
}

namespace OpenHurricane {
    template <> void geometryArray<symmTensor, cellMesh>::writeOutput(fileOsstream &fos) const {
        if (argParse::checkWriteSolution()) {
            bool isok = true;

            for (integer n = 0; n < this->internalArraySize(); ++n) {
                if (isnan(this->operator[](n).xx()) || isinf(this->operator[](n).xx()) ||
                    isnan(this->operator[](n).xy()) || isinf(this->operator[](n).xy()) ||
                    isnan(this->operator[](n).xz()) || isinf(this->operator[](n).xz()) ||
                    isnan(this->operator[](n).yy()) || isinf(this->operator[](n).yy()) ||
                    isnan(this->operator[](n).yz()) || isinf(this->operator[](n).yz()) ||
                    isnan(this->operator[](n).zz()) || isinf(this->operator[](n).zz())) {
                    isok = false;
                    break;
                }
            }
            if (!isok) {
                LFatal("Not a number or infinity of %s", this->name().c_str());
            }
        }

        // Only write the internale filed's value.
        Array<symmTensor>::writeToStream(fos, this->internalArraySize());
    }

    template <>
    void
    OpenHurricane::geometryArray<symmTensor, cellMesh>::writeOutputByMaster(fileOsstream &fos) const {
        if (HurMPI::parRun()) {
            integerList nSizeL;
            integerList displs;
            integer allSize = 0;
            realArray rootF;
            realArray cmpF;

            nSizeL.resize(HurMPI::getProcSize(), Zero);
            nSizeL[HurMPI::getProcRank()] = this->internalArraySize();
            if (HurMPI::master()) {
                displs.resize(HurMPI::getProcSize(), Zero);
            }
            HurMPI::gatherList(nSizeL, HurMPI::masterNo(), HurMPI::getComm());

            if (HurMPI::master()) {
                for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                    displs[ip] = displs[ip - 1] + nSizeL[ip - 1];
                }
                for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                    allSize += nSizeL[ip];
                }
            }
            cmpF.resize(this->internalArraySize());
            if (HurMPI::master()) {
                rootF.resize(allSize);
            }
            if (argParse::checkWriteSolution()) {
                bool isok = true;

                for (integer n = 0; n < this->internalArraySize(); ++n) {
                    if (isnan(this->operator[](n).xx()) || isinf(this->operator[](n).xx()) ||
                        isnan(this->operator[](n).xy()) || isinf(this->operator[](n).xy()) ||
                        isnan(this->operator[](n).xz()) || isinf(this->operator[](n).xz()) ||
                        isnan(this->operator[](n).yy()) || isinf(this->operator[](n).yy()) ||
                        isnan(this->operator[](n).yz()) || isinf(this->operator[](n).yz()) ||
                        isnan(this->operator[](n).zz()) || isinf(this->operator[](n).zz())) {
                        isok = false;
                        break;
                    }
                }
                if (!isok) {
                    LFatal("Not a number or infinity of %s", this->name().c_str());
                }
            }

            HurMPI::barrier(HurMPI::getComm());

            for (integer i = 0; i < feature<symmTensor>::nElements_; ++i) {
                for (integer j = 0; j < this->internalArraySize(); ++j) {
                    cmpF[j] = this->operator[](j)[i];
                }

                HurMPI::Request request;
                HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                                 rootF.data(), nSizeL.data(), displs.data(), feature<real>::MPIType,
                                 HurMPI::masterNo(), HurMPI::getComm(), &request);
                HurMPI::wait(&request, MPI_STATUSES_IGNORE);

                if (HurMPI::master()) {
                    fos.write(reinterpret_cast<const char *>(&rootF[0]), rootF.byteSize());
                }
            }
        } else {
            writeOutput(fos);
        }
    }

    template <>
    void OpenHurricane::geometryArray<symmTensor, cellMesh>::writeOutput(fileOsstream &fos,
                                                                     const integer fzid) const {
        const auto &displs = mesh().globalFaceZoneInfo(fzid).faceDispls();
        const auto &recvcnt = mesh().globalFaceZoneInfo(fzid).faceRecvcnt();
        const auto &faces = mesh().faces();
        const auto &fw = mesh().faceWgt();
        const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();
        const auto fsize = mesh().globalFaceZoneInfo(fzid).fZ().size();
        const auto nTotalFaces = mesh().globalFaceZoneInfo(fzid).totalFaces();

        realArray faceValue(fsize, Zero);
        realArray writeField;
        if (HurMPI::parRun()) {
            if (HurMPI::master()) {
                writeField.resize(nTotalFaces, Zero);
            }
        }
        for (int j = 0; j < symmTensor::nElements_; ++j) {
            integer count = 0;
            for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
                const integer cl = faces[fi].leftCell();
                const integer cr = faces[fi].rightCell();
                faceValue[count++] =
                    fw[fi] * this->operator[](cl)[j] + (1.0 - fw[fi]) * this->operator[](cr)[j];
            }
            if (!HurMPI::parRun()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < fsize) {
                        sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                             << " ";
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        if (i < fsize) {
                            sstr << std::setprecision(feature<real>::precision) << faceValue[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (fsize) {
                        fos.write(reinterpret_cast<const char *>(&faceValue[0]),
                                  fsize * sizeof(real));
                    }
                }
                continue;
            }

            HurMPI::Request request;
            HurMPI::igatherv(faceValue.data(), fsize, feature<real>::MPIType, writeField.data(),
                             recvcnt.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                if (fos.format() == IOsstream::ASCII_FORMAT) {
                    std::stringstream sstr;
                    integer i = 0;
                    while (i < writeField.size()) {
                        sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                             << " ";
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        if (i < writeField.size()) {
                            sstr << std::setprecision(feature<real>::precision) << writeField[i++]
                                 << " ";
                        }
                        sstr << "\n";
                    }
                    fos.os() << sstr.str().c_str();
                } else {
                    if (writeField.size()) {
                        fos.write(reinterpret_cast<const char *>(&writeField[0]),
                                  writeField.size() * sizeof(real));
                    }
                }
            }
        }
    }

} // namespace OpenHurricane

template <>
void OpenHurricane::geometryArray<OpenHurricane::symmTensor, OpenHurricane::cellMesh>::writeMinMaxOutput(
    fileOsstream &fos, const integer fzid) const {
    const auto &faces = mesh().faces();
    const auto &fw = mesh().faceWgt();
    const auto &fZ = mesh().globalFaceZoneInfo(fzid).fZ();

    for (integer j = 0; j < symmTensor::nElements_; ++j) {
        real minX = veryLarge;
        real maxX = -veryLarge;
        for (integer fi = fZ.firstIndex(); fi <= fZ.lastIndex(); ++fi) {
            const integer cl = faces[fi].leftCell();
            const integer cr = faces[fi].rightCell();
            real vv = 0.5 * this->operator[](cl)[j] + 0.5 * this->operator[](cr)[j];

            minX = min(minX, vv);
            maxX = max(maxX, vv);
        }
        if (HurMPI::parRun()) {
            HurMPI::reduce(minX, MPI_MIN);
            HurMPI::reduce(maxX, MPI_MAX);
        }
        if (HurMPI::master()) {
            fos.write(reinterpret_cast<const char *>(&minX), sizeof(double));
            fos.write(reinterpret_cast<const char *>(&maxX), sizeof(double));
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::symmTensor, OpenHurricane::cellMesh>::writeRelay(
    hdf5O &fos, const bool writeLast, const bool writeToGroup) const {
    integerList nSizeL;
    integerList displs;
    integer allSize;
    symmTensorArray rootV;
    realArray rootF;
    realArray cmpF;
    if (HurMPI::parRun()) {
        nSizeL.resize(HurMPI::getProcSize(), Zero);
        nSizeL[HurMPI::getProcRank()] = this->internalArraySize();
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
        cmpF.resize(this->internalArraySize());
        if (HurMPI::master()) {
            rootF.resize(allSize);
            rootV.resize(allSize);
        }

        HurMPI::barrier(HurMPI::getComm());

        for (integer i = 0; i < feature<symmTensor>::nElements_; ++i) {
            for (integer j = 0; j < this->internalArraySize(); ++j) {
                cmpF[j] = this->operator[](j)[i];
            }

            HurMPI::Request request;
            HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                             rootF.data(), nSizeL.data(), displs.data(), feature<real>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            for (integer j = 0; j < rootF.size(); ++j) {
                rootV[j][i] = rootF[j];
            }
        }
        if (HurMPI::master()) {
            if (writeLast || writeToGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.write(rootV, real(1), rootV.size(), groupName, this->name());
                for (integer i = 0; i < feature<symmTensor>::nElements_; ++i) {
                    string rhsName = "rhs0" + toString(i);
                    fos.writeRealAttributeToDataset(rhs0_[i], rhsName, groupName, this->name());
                }
            } else {
                fos.write(rootV, real(1), rootV.size(), this->name());
                for (integer i = 0; i < feature<symmTensor>::nElements_; ++i) {
                    string rhsName = "rhs0" + toString(i);
                    fos.writeRealAttributeToDataset(rhs0_[i], rhsName, this->name());
                }
            }
        }
    } else {
        if (writeLast || writeToGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.write(*this, real(1), this->internalArraySize(), groupName, this->name());

            for (integer i = 0; i < feature<symmTensor>::nElements_; ++i) {
                string rhsName = "rhs0" + toString(i);
                fos.writeRealAttributeToDataset(rhs0_[i], rhsName, groupName, this->name());
            }
        } else {
            fos.write(*this, real(1), this->internalArraySize(), this->name());

            for (integer i = 0; i < feature<symmTensor>::nElements_; ++i) {
                string rhsName = "rhs0" + toString(i);
                fos.writeRealAttributeToDataset(rhs0_[i], rhsName, this->name());
            }
        }
    }
    if (writeLast) {
        if (lastArrayArray_.size() == 0) {
            return;
        }
        if (HurMPI::parRun()) {
            for (integer k = 0; k < lastArrayArray_.size(); ++k) {
                HurMPI::barrier(HurMPI::getComm());
                for (integer i = 0; i < feature<symmTensor>::nElements_; ++i) {
                    for (integer j = 0; j < this->internalArraySize(); ++j) {
                        cmpF[j] = lastArrayArray_[k][j][i];
                    }

                    HurMPI::Request request;
                    HurMPI::igatherv(cmpF.data(), this->internalArraySize(), feature<real>::MPIType,
                                     rootF.data(), nSizeL.data(), displs.data(),
                                     feature<real>::MPIType, HurMPI::masterNo(), HurMPI::getComm(),
                                     &request);
                    HurMPI::wait(&request, MPI_STATUSES_IGNORE);

                    for (integer j = 0; j < rootF.size(); ++j) {
                        rootV[j][i] = rootF[j];
                    }
                }
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(k + 1);
                    fos.write(rootV, rootV.size(), groupName, this->name());
                }
            }
        } else {
            for (integer k = 0; k < lastArrayArray_.size(); ++k) {
                string groupName = "timeGroupm";
                groupName += toString(k + 1);
                fos.write(lastArrayArray_[k], this->internalArraySize(), groupName, this->name());
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::symmTensor, OpenHurricane::cellMesh>::readRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {
    if (HurMPI::parRun()) {
        symmTensorArray allTmpData;
        symmTensorArray relayData;
        if (HurMPI::master()) {
            if (readLast || readFromGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.read(allTmpData, groupName, this->name());
                for (integer i = 0; i < feature<symmTensor>::nElements_; i++) {
                    string rhsName = "rhs0" + toString(i);
                    fos.readRealAttributeFromDataset(rhs0_[i], rhsName, groupName, this->name());
                }
            } else {
                fos.read(allTmpData, this->name());
                for (integer i = 0; i < feature<symmTensor>::nElements_; i++) {
                    string rhsName = "rhs0" + toString(i);
                    fos.readRealAttributeFromDataset(rhs0_[i], rhsName, this->name());
                }
            }
            const std::map<integer, integer> &indexMap = mesh().indexMap();
            const auto &originCellIndex = mesh().originCellIndex();
            relayData.resize(allTmpData.size());

            integer offset = 0;

            for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                    std::map<integer, integer>::const_iterator iter =
                        indexMap.find(originCellIndex[czi][i]);
                    relayData[i + offset] = allTmpData[iter->second];
                }
                offset += originCellIndex[czi].size();
            }
        }
        symmTensorArray tmpData;
        relayScatterFunc::relayScatter<symmTensor>(mesh().cellOfProc(), relayData, tmpData);
        relayScatterFunc::relayReorder<symmTensor>(tmpData, *this, mesh().perm(), mesh().iperm());
        integer ng = 0;
        if (HurMPI::master()) {
            if (!fos.exist("nTimeGroups")) {
                fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
            }
        }
        HurMPI::bcast(&ng, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(i + 1);
                    fos.read(allTmpData, groupName, this->name());
                    const std::map<integer, integer> &indexMap = mesh().indexMap();
                    const auto &originCellIndex = mesh().originCellIndex();

                    integer offset = 0;

                    for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                        for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                            std::map<integer, integer>::const_iterator iter =
                                indexMap.find(originCellIndex[czi][i]);
                            relayData[i + offset] = allTmpData[iter->second];
                        }
                        offset += originCellIndex[czi].size();
                    }
                }
                relayScatterFunc::relayScatter<symmTensor>(mesh().cellOfProc(), relayData, tmpData);
                relayScatterFunc::relayReorder<symmTensor>(tmpData, lastArrayArray_[i],
                                                           mesh().perm(), mesh().iperm());
            }
        }
    } else {
        symmTensorArray tmpData;
        symmTensorArray relayData;
        if (readLast || readFromGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.read(tmpData, groupName, this->name());
            for (integer i = 0; i < feature<symmTensor>::nElements_; i++) {
                string rhsName = "rhs0" + toString(i);
                fos.readRealAttributeFromDataset(rhs0_[i], rhsName, groupName, this->name());
            }
        } else {
            fos.read(tmpData, this->name());
            for (integer i = 0; i < feature<symmTensor>::nElements_; i++) {
                string rhsName = "rhs0" + toString(i);
                fos.readRealAttributeFromDataset(rhs0_[i], rhsName, this->name());
            }
        }
        const std::map<integer, integer> &indexMap = mesh().indexMap();
        const auto &originCellIndex = mesh().originCellIndex();
        relayData.resize(tmpData.size());

        integer offset = 0;

        for (integer czi = 0; czi < originCellIndex.size(); czi++) {
            for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                std::map<integer, integer>::const_iterator iter =
                    indexMap.find(originCellIndex[czi][i]);
                relayData[i + offset] = tmpData[iter->second];
            }
            offset += originCellIndex[czi].size();
        }
        relayScatterFunc::relayReorder<symmTensor>(relayData, *this, mesh().perm(), mesh().iperm());
        integer ng = 0;
        if (!fos.exist("nTimeGroups")) {
            fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
        }
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                symmTensorArray tmpData;
                string groupName = "timeGroupm";
                groupName += toString(i + 1);
                fos.read(tmpData, groupName, this->name());
                const std::map<integer, integer> &indexMap = mesh().indexMap();
                const auto &originCellIndex = mesh().originCellIndex();

                integer offset = 0;

                for (integer czi = 0; czi < originCellIndex.size(); czi++) {
                    for (integer i = 0; i < originCellIndex[czi].size(); i++) {
                        std::map<integer, integer>::const_iterator iter =
                            indexMap.find(originCellIndex[czi][i]);
                        relayData[i + offset] = tmpData[iter->second];
                    }
                    offset += originCellIndex[czi].size();
                }
                relayScatterFunc::relayReorder<symmTensor>(relayData, lastArrayArray_[i],
                                                           mesh().perm(), mesh().iperm());
            }
        }
    }
}

template <>
void OpenHurricane::geometryArray<OpenHurricane::symmTensor, OpenHurricane::cellMesh>::interpolateRelay(
    const hdf5I &fos, const bool readLast, const bool readFromGroup) {
    if (HurMPI::parRun()) {
        //symmTensorArray allTmpData;
        symmTensorArray relayData;
        if (HurMPI::master()) {
            if (readLast || readFromGroup) {
                string groupName = "timeGroup";
                groupName += toString(0);
                fos.read(relayData, groupName, this->name());
                for (integer i = 0; i < feature<symmTensor>::nElements_; i++) {
                    string rhsName = "rhs0" + toString(i);
                    fos.readRealAttributeFromDataset(rhs0_[i], rhsName, groupName, this->name());
                }
            } else {
                fos.read(relayData, this->name());
                for (integer i = 0; i < feature<symmTensor>::nElements_; i++) {
                    string rhsName = "rhs0" + toString(i);
                    fos.readRealAttributeFromDataset(rhs0_[i], rhsName, this->name());
                }
            }
        }
        symmTensorArray tmpData;
        relayScatterFunc::relayScatter<symmTensor>(mesh().tarOfProc(), relayData, tmpData);

        const auto &nbr = mesh().sorKNN();
        const auto &sor = mesh().sorCellCentre();
        const auto &tar = mesh().cellCentre();

        meshInterpolation itp(nbr, sor, tar);
        itp.interpolate<symmTensor>(tmpData, *this);

        integer ng = 0;
        if (HurMPI::master()) {
            if (!fos.exist("nTimeGroups")) {
                fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
            }
        }
        HurMPI::bcast(&ng, 1, feature<integer>::MPIType, HurMPI::masterNo(), HurMPI::getComm());
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                if (HurMPI::master()) {
                    string groupName = "timeGroupm";
                    groupName += toString(i + 1);
                    fos.read(relayData, groupName, this->name());
                }
                relayScatterFunc::relayScatter<symmTensor>(mesh().tarOfProc(), relayData, tmpData);
                itp.interpolate<symmTensor>(tmpData, lastArrayArray_[i]);
            }
        }
    } else {
        symmTensorArray tmpData;
        symmTensorArray relayData;
        if (readLast || readFromGroup) {
            string groupName = "timeGroup";
            groupName += toString(0);
            fos.read(tmpData, groupName, this->name());
            for (integer i = 0; i < feature<symmTensor>::nElements_; i++) {
                string rhsName = "rhs0" + toString(i);
                fos.readRealAttributeFromDataset(rhs0_[i], rhsName, groupName, this->name());
            }
        } else {
            fos.read(tmpData, this->name());
            for (integer i = 0; i < feature<symmTensor>::nElements_; i++) {
                string rhsName = "rhs0" + toString(i);
                fos.readRealAttributeFromDataset(rhs0_[i], rhsName, this->name());
            }
        }

        relayScatterFunc::relayReorder<symmTensor>(relayData, *this, mesh().perm(), mesh().iperm());

        const auto &nbr = mesh().sorKNN();
        const auto &sor = mesh().sorCellCentre();
        const auto &tar = mesh().cellCentre();

        meshInterpolation itp(nbr, sor, tar);
        itp.interpolate<symmTensor>(tmpData, *this);

        integer ng = 0;
        if (!fos.exist("nTimeGroups")) {
            fos.readIntegerAttributeFromFile(ng, "nTimeGroups");
        }
        if (readLast && ng > 0) {
            const auto minS = min(lastArrayArray_.size(), ng);
            for (integer i = 0; i < minS; ++i) {
                string groupName = "timeGroupm";
                groupName += toString(i + 1);
                fos.read(relayData, groupName, this->name());

                relayScatterFunc::relayScatter<symmTensor>(mesh().tarOfProc(), relayData, tmpData);
                itp.interpolate<symmTensor>(tmpData, lastArrayArray_[i]);
            }
        }
    }
}