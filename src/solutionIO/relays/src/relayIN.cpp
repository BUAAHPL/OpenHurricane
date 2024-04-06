/*!
 * \file relayIN.cpp
 * \brief The subroutines and functions of reading relay files
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
#include "relayIN.hpp"
#include "controllerSwitch.hpp"
#include "geometryMesh.hpp"
#include "hdf5I.hpp"
#include "hdf5O.hpp"
#include "iteration.hpp"

OpenHurricane::relayIN::relayIN(iteration &iter)
    : relay(), iter_(iter), isInterpolation_(false), steadyOrUnsteady_(false) {
    if (iter.cont().found("iteration")) {
        const auto &iterCont = iter.cont().subController("iteration");
        if (iterCont.found("relayType")) {
            isInterpolation_ = controllerSwitch(iterCont)("relayType", isInterpolation_,
                                                          "interpolationData", "originData");
        }
    }
}

void OpenHurricane::relayIN::readRelay(const fileName &restartFrom) {
    if (!isInterpolation_) {
        PLInfo("\n    Reading relay file from: %s\n\n", restartFrom.c_str());
    } else {
        PLInfo("\n    Reading and interpolating from relay file: %s\n\n", restartFrom.c_str());
    }
    hdf5I myh5(restartFrom);
    myh5.open();
    readRelay(myh5, restartFrom);
    myh5.close();
}

void OpenHurricane::relayIN::readRelay(const hdf5I &fos, const fileName &restartFrom) {
    const geometryMesh &mesh = iter_.findObject<geometryMesh>(iter_.name());
    auto nv = getIntegerAttributeFromFile(fos, attN::nVariables);
    auto totalStep = getIntegerAttributeFromFile(fos, attN::totalStep);
    auto state = getIntegerAttributeFromFile(fos, attN::state);

    auto progName = getStringAttributeFromFile(fos, attN::program);
    auto fileVersion = getStringAttributeFromFile(fos, attN::version);
    auto fileDataTime = getStringAttributeFromFile(fos, attN::dateTime);
    auto meshFile = getStringAttributeFromFile(fos, attN::meshFile);

    auto nTimeGroups = getIntegerAttributeFromFile(fos, attN::nTimeGroups);

    iter_.setTotalStep(totalStep);
    iter_.setNTimeGroups(nTimeGroups);

    if (state != 0) {
        getPhysicalTimeInfoFromFile(fos);
    }
    stringList nvN = getVarName(fos);

    bool readFromGroup = false;
    integer ng = getIntegerAttributeFromFile(fos, attN::nTimeGroups);

    if (state != 0) {
        readFromGroup = true;
    }
    bool readLast = false;
    if (state != 0) {
        steadyOrUnsteady_ = true;
        readLast = iter_.isReadLastFromRelay() && (ng > 0);
    }
    if (!isInterpolation_) {
        createIndexMap(fos);
        const_cast<geometryMesh &>(mesh).readRelay(fos, nvN, readLast, readFromGroup);
    } else {
        if (const_cast<geometryMesh &>(mesh).getInterpolationSource(fos)) {
            const_cast<geometryMesh &>(mesh).interpolateRelay(fos, nvN, readLast, readFromGroup);
        } else {
            LFatal("Cannot found cell center info of file: %s for interpolating new solution. "
                   "Please check!",
                   restartFrom.c_str());
        }
    }
}

void OpenHurricane::relayIN::createIndexMap(const hdf5I &fos) {
    geometryMesh &mesh = iter_.findObjectRef<geometryMesh>(iter_.name());
    integer nCells = mesh.nCells();
    HurMPI::reduce(nCells, MPI_SUM);

    if (HurMPI::master()) {
        integerList originIndex(mesh.nCells());
        std::map<integer, integer> &indexMap = mesh.indexMap();

        fos.read(originIndex, "cellOriginIndex");

        for (integer i = 0; i < nCells; i++) {
            indexMap[originIndex[i]] = i;
        }
    }
}

hur_nodiscard OpenHurricane::integer
OpenHurricane::relayIN::getIntegerAttributeFromFile(const hdf5I &fos,
                                                    const string &attrName) const {
    integer iatrr = 0;
    if (HurMPI::master()) {
        fos.readIntegerAttributeFromFile(iatrr, attrName);
    }
    HurMPI::bcast(&iatrr, 1, feature<integer>::MPIType);

    return iatrr;
}

hur_nodiscard std::string
OpenHurricane::relayIN::getStringAttributeFromFile(const hdf5I &fos, const string &attrName) const {
    std::string mstr;
    if (HurMPI::master()) {
        fos.readStringAttributeFromFile(mstr, attrName);
    }
    HurMPI::bcastString(mstr);
    return mstr;
}

hur_nodiscard OpenHurricane::real
OpenHurricane::relayIN::getRealAttributeFromFile(const hdf5I &fos, const string &attrName) const {
    real ratrr = 0;
    if (HurMPI::master()) {
        fos.readRealAttributeFromFile(ratrr, attrName);
    }
    HurMPI::bcast(&ratrr, 1, feature<real>::MPIType);

    return ratrr;
}

hur_nodiscard bool OpenHurricane::relayIN::fileExist(const hdf5I &fos,
                                                     const string &attrName) const {
    bool lastTimeStepExist = false;
    if (HurMPI::master()) {
        if (fos.exist(attrName)) {
            lastTimeStepExist = true;
        }
    }
    HurMPI::bcast(&lastTimeStepExist, 1, feature<bool>::MPIType);

    return lastTimeStepExist;
}

hur_nodiscard OpenHurricane::realArray
OpenHurricane::relayIN::getRealArrayFromFile(const hdf5I &fos, const string &attrName) const {
    realArray rarray;
    if (HurMPI::master()) {
        fos.read(rarray, attrName);
    }
    integer si = rarray.size();
    HurMPI::bcast(&si, 1, feature<integer>::MPIType);
    if (!HurMPI::master()) {
        rarray.resize(si, Zero);
    }
    HurMPI::bcastList(rarray);

    return rarray;
}

hur_nodiscard std::string OpenHurricane::relayIN::getStringFromFile(const hdf5I &fos,
                                                                    const string &attrName) const {
    std::string nameStr;
    if (HurMPI::master()) {
        fos.readString(nameStr, attrName);
    }
    HurMPI::bcastString(nameStr);
    return nameStr;
}

void OpenHurricane::relayIN::setPhysicalTimeInfo(const real totalPhyTimes,
                                                 const realArray &lastTimeStep) {
    if (iter_.hasPhysicalTimeStep()) {
        iter_.pTStep().setTotalTime(totalPhyTimes);
        if (lastTimeStep.size() != 0) {
            if (iter_.pTStep().lastTimeStep().size() < lastTimeStep.size()) {
                iter_.pTStep().lastTimeStep().resize(lastTimeStep.size());
            }
            for (integer i = 0; i < lastTimeStep.size(); ++i) {
                iter_.pTStep().lastTimeStep()[i] = lastTimeStep[i];
            }
        }
    }
}

void OpenHurricane::relayIN::getPhysicalTimeInfoFromFile(const hdf5I &fos) {
    auto ptimes = getRealAttributeFromFile(fos, attN::time);
    realArray lastTimeStep;
    auto lastTimeStepExist = fileExist(fos, attN::lastTimeStep);
    if (lastTimeStepExist) {
        lastTimeStep = getRealArrayFromFile(fos, attN::lastTimeStep);
    }

    setPhysicalTimeInfo(ptimes, lastTimeStep);
}

hur_nodiscard OpenHurricane::stringList OpenHurricane::relayIN::getVarName(const hdf5I &fos) const {
    stringList nvN;
    std::string nameStr = getStringFromFile(fos, attN::variableName);
    split(nameStr, nvN, ",");
    return nvN;
}
