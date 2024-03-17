/*!
 * \file relayOUT.cpp
 * \brief The subroutines and functions of writting relay files
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
#include "relayOUT.hpp"
#include "casePredifine.hpp"
#include "geometryMesh.hpp"
#include "hdf5I.hpp"
#include "hdf5O.hpp"
#include "iteration.hpp"

OpenHurricane::relayOUT::relayOUT(const iteration &iter) : relay(), iter_(iter) {}

void OpenHurricane::relayOUT::writeRelay(const fileName &relayFileName,
                                         const fileName &meshFileName) const {
    const geometryMesh &mesh = iter_.findObject<geometryMesh>(iter_.name());

    if (HurMPI::master()) {
        ofos_.open(relayFileName);
    }
    if (HurMPI::master()) {
        writeStringAttriToFile(programName::name, attN::program);
        writeStringAttriToFile(programName::getVersion(), attN::version);
        writeStringAttriToFile(clock::dateTime(), attN::dateTime);
        writeStringAttriToFile(meshFileName.name(), attN::meshFile);

        writeIntAttriToFile(HurMPI::getProcSize(), attN::nProcessor);

        // 0 - only data
        // 1 - with grid and data
        // 2 - only grid
        // 3 - case config and grid
        writeIntAttriToFile(relayDefine::ONLY_DATA, attN::dataType);
        writeIntAttriToFile(iter_.totalStep(), attN::totalStep);
        if (iter_.isSteadyFlow()) {
            ofos_.writeRealAttributeToFile(Zero, attN::time);
            writeIntAttriToFile(0, attN::state);
        } else {
            ofos_.writeRealAttributeToFile(iter_.pTStep().totalTime(), attN::time);
            writeIntAttriToFile(1, attN::state);
            ofos_.write(iter_.pTStep().lastTimeStep(), attN::lastTimeStep);
        }
        const auto nVar = mesh.outputRelaySize();
        const auto varName = mesh.outputRelayNameDocList();

        writeIntAttriToFile(nVar, attN::nVariables);
        ofos_.write(varName, attN::variableName);
    }
    /*string gridGroupName = "grid";
    mesh.writeMesh(myh5, gridGroupName);*/

    writeOriIndex(ofos_);
    writeCellCentre(ofos_);

    if (iter_.isSteadyFlow()) {
        mesh.writeRelay(ofos_, false, false);
    } else {
        mesh.writeRelay(ofos_, iter_.isWriteLastToRelay(), true);
    }

    ofos_.close();
}

void OpenHurricane::relayOUT::writeOriIndex(hdf5O &fos) const {
    integerList nSizeL;
    integerList displs;
    integer allSize = 0;
    const geometryMesh &mesh = iter_.findObject<geometryMesh>(iter_.name());
    if (HurMPI::parRun()) {
        nSizeL.resize(HurMPI::getProcSize(), Zero);
        nSizeL[HurMPI::getProcRank()] = mesh.nCells();
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

        integerArray tmpData(mesh.nCells());
        for (integer i = 0; i < mesh.nCells(); i++) {
            tmpData[i] = mesh.cells()[i].orignIndex();
        }

        integerArray rootF;
        if (HurMPI::master()) {
            rootF.resize(allSize);
        }

        HurMPI::barrier(HurMPI::getComm());

        HurMPI::Request request;
        HurMPI::igatherv(tmpData.data(), mesh.nCells(), feature<integer>::MPIType, rootF.data(),
                         nSizeL.data(), displs.data(), feature<integer>::MPIType,
                         HurMPI::masterNo(), HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);

        if (HurMPI::master()) {
            fos.write(rootF, attN::cellOriginIndex);
        }
    } else {
        integerArray expData(mesh.nCells());
        for (integer i = 0; i < mesh.nCells(); i++) {
            expData[i] = mesh.cells()[i].orignIndex();
        }
        fos.write(expData, attN::cellOriginIndex);
    }
}

void OpenHurricane::relayOUT::writeCellCentre(hdf5O &fos) const {
    const geometryMesh &mesh = iter_.findObject<geometryMesh>(iter_.name());

    auto cellCentre = mesh.allCellCentre();
    if (HurMPI::master()) {
        fos.write(cellCentre, attN::cellCentre);
    }
}
