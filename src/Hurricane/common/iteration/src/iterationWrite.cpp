/*!
 * \file iterationWrite.cpp
 * \brief The subroutines and functions of CFD time advance iteration
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
#include "casePredifine.hpp"
#include "geometryMesh.hpp"
#include "hdf5I.hpp"
#include "hdf5O.hpp"
#include "iteration.hpp"
#include "relayOUT.hpp"
#include "solutionWrite.hpp"
#include "tecplotWriter.hpp"
#include "writeFaceZone.hpp"

void OpenHurricane::iteration::writeDualTimePhyInfo() const {
    if (HurMPI::master()) {
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::left << std::setfill(' ') << std::setw(10) << "Iter";
        std::cout << std::left << std::setfill(' ') << std::setw(20) << "Current Time(s)";
        std::cout << std::left << std::setfill(' ') << std::setw(20) << "Total Time(s)";
        std::cout << std::left << std::setfill(' ') << std::setw(20) << "Time step(s)";
        std::cout << std::left << std::setfill(' ') << std::setw(20) << "Remaining" << std::endl;
        std::cout << std::left << std::setfill(' ') << std::setw(10) << totalStep();
        std::cout.setf(std::ios::showpoint);
        std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(5)
                  << pTStep().totalTime();
        std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(5)
                  << pTStep().maxTime();
        std::cout << std::left << std::setfill(' ') << std::setw(20) << std::setprecision(5)
                  << pTStep().pTimeStep();
        std::cout.unsetf(std::ios::showpoint);
        std::cout << std::left << std::setfill(' ') << std::setw(20) << maxStep() - cStep()
                  << std::endl;
        std::cout << std::endl;
        std::cout << std::right;
    }
}

void OpenHurricane::iteration::write() const {
    monitoring();

    if (writeFaceZonePtr_) {
        writeFaceZonePtr_->write();
    }

    if (solWritePtr_) {
        solWritePtr_->write();
    }

    if (cStep_ % writeOutputStep_ == 0 || (end() && !argParse::noWriteAtEnd())) {
        if (isGUIMode()) {
            const auto caseN = getCaseConfigName();
            const auto outN = getRelayFileName();
            writeCaseConfig(caseN);
            writeRelay(outN, caseN);
        } else {
            writeRelay(getRelayFileName(), meshFileName());
        }

        if (cont().subController("iteration").found("cellLoadWeight")) {
            const geometryMesh &mesh = this->findObject<geometryMesh>(name());
            if (mesh.hasCellLoadWeights()) {
                fileName outN = meshFileName();
                auto pathOut = outN.parentPath();
                auto fname = outN.name(true);
                string fext;
                fext += outN.ext();

                fname += "-";
                fname += toString(totalStep_);
                outN = fname + fext;
                outN = pathOut / outN;

                writeMesh(outN);
            } else {
                PLWarning("Cannot write cell-load-weights");
            }
        }
    }
}

void OpenHurricane::iteration::writeMesh(const fileName &meshFileName) const {
    const geometryMesh &mesh = this->findObject<geometryMesh>(name());
    hdf5O myh5(meshFileName);
    if (HurMPI::master()) {
        myh5.open();

        myh5.writeStringAttributeToFile(programName::name, string("program"));
        myh5.writeStringAttributeToFile(programName::getVersion(), string("version"));
        myh5.writeStringAttributeToFile(clock::dateTime(), string("dateTime"));
        myh5.writeStringAttributeToFile(meshFileName.name(), string("meshFile"));
        myh5.writeIntegerAttributeToFile(HurMPI::getProcSize(), string("nProcessor"));

        // 0 - only data
        // 1 - with grid and data
        // 2 - only grid
        // 3 - grid and case config
        myh5.writeIntegerAttributeToFile(Options::writeRelayOption::ONLY_GRID, string("dataType"));
    }
    string gridGroupName = "grid";
    mesh.writeMesh(myh5, gridGroupName);
    myh5.close();
}

void OpenHurricane::iteration::writeCaseConfig(const fileName &meshFileName) const {
    const geometryMesh &mesh = this->findObject<geometryMesh>(name());
    hdf5O myh5(meshFileName);
    if (HurMPI::master()) {
        myh5.open();

        myh5.writeStringAttributeToFile(programName::name, string("program"));
        myh5.writeStringAttributeToFile(programName::getVersion(), string("version"));
        myh5.writeStringAttributeToFile(clock::dateTime(), string("dateTime"));
        myh5.writeStringAttributeToFile(meshFileName.name(), string("meshFile"));
        myh5.writeIntegerAttributeToFile(HurMPI::getProcSize(), string("nProcessor"));

        // 0 - only data
        // 1 - with grid and data
        // 2 - only grid
        // 3 - grid and case config
        myh5.writeIntegerAttributeToFile(Options::writeRelayOption::CASE_CONFIG,
                                         string("dataType"));
    }
    string gridGroupName = "grid";
    mesh.writeMesh(myh5, gridGroupName);
    if (HurMPI::master()) {
        string controllerGroupName = "controller";
        std::stringstream controllerString;
        cont().write(controllerString);
        myh5.writeString(controllerString.str(), controllerGroupName);
    }
    myh5.close();
}

void OpenHurricane::iteration::writeRelay(const fileName &relayFileName,
                                      const fileName &meshFileName) const {
    relayOUT myOut(*this);
    myOut.writeRelay(relayFileName, meshFileName);
}

void OpenHurricane::iteration::writeResult() const {}

void OpenHurricane::iteration::writeFieldVarResult(const fileName &outN,
                                               const stringList &varList) const {
    Pout << "    Writting result to file: " << outN << std::endl;
    const geometryMesh &mesh = this->findObject<geometryMesh>(name());
    if (argParse::isASCII()) {
        fileOsstream fos(outN);
        mesh.writeOutputToTecplot(fos, varList);
        fos.close();
    } else {
        fileOsstream fos(outN, IOsstream::BINARY_FORMAT);
        if (HurMPI::multiNodes() || argParse::tecplotFileWriteOnMaster()) {
            tecplotWriter::writeResultToTecplotByMaster(fos, mesh, varList);
        } else {
            tecplotWriter::writeResultToTecplot(fos, mesh, varList);
        }
        fos.close();
    }
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::iteration::setFileName(const string &type) const {
    fileName outN = outputName();
    auto pathOut = outN.parentPath();
    auto fname = outN.name(true);
    string fext;
    if (argParse::isASCII()) {
        fext += outN.ext();
    } else {
        fext = ".plt";
    }
    fname += "-";
    fname += type;
    fname += "-";
    fname += toString(totalStep_);
    outN = fname + fext;
    outN = pathOut / outN;
    return outN;
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::iteration::getRelayFileName() const {
    fileName outN = dataName();
    auto pathOut = outN.parentPath();
    auto fname = outN.name(true);
    string fext;
    fext += outN.ext();

    fname += "-";
    fname += toString(totalStep_);
    outN = fname + fext;
    outN = pathOut / outN;
    return outN;
}

hur_nodiscard OpenHurricane::fileName OpenHurricane::iteration::getCaseConfigName() const {
    fileName outN = caseFileName();
    auto pathOut = outN.parentPath();
    auto fname = outN.name(true);
    string fext;
    fext += outN.ext();

    fname += "-";
    fname += toString(totalStep_);
    outN = fname + fext;
    outN = pathOut / outN;
    return outN;
}

void OpenHurricane::iteration::writeSubIterResiduals() const {
    if (!hasSubIteration()) {
        return;
    }
    subMonitoring();
}

void OpenHurricane::iteration::writeMeshToTecplot() const {
    const geometryMesh &mesh = this->findObject<geometryMesh>(name());

    fileName meshFile = caseFileName().name(true);
    
    meshFile += "_Grid.dat";
    if (argParse::isASCII()) {
        fileOsstream fos(meshFile);
        mesh.writeMeshToTecplot(fos);
        fos.close();
    } else {
        fileOsstream fos(meshFile, IOsstream::BINARY_FORMAT);
        tecplotWriter::writeMesh(fos, mesh);
        fos.close();
    }
}

void OpenHurricane::iteration::setSolWrite(const flowModel &flows) {
    solWritePtr_.clear();
    const auto &wCont = cont().subController("writeControl");
    const auto &swCont = wCont.subController("solutionWrite");
    solWritePtr_ = solutionWrite::creator(flows, *this, flows.mesh(), swCont);
}

void OpenHurricane::iteration::setWriteFaceZonePtr(writeFaceZone *newPtr) {
    writeFaceZonePtr_.reset(newPtr);
    newPtr = nullptr;
}

hur_nodiscard OpenHurricane::solutionWrite &OpenHurricane::iteration::solWrite() noexcept {
    if (!solWritePtr_) {
        LFatal("Attempt to access null pointer");
    }
    return *solWritePtr_;
}

hur_nodiscard const OpenHurricane::solutionWrite &OpenHurricane::iteration::solWrite() const noexcept {
    if (!solWritePtr_) {
        LFatal("Attempt to access null pointer");
    }
    return *solWritePtr_;
}
