/*!
 * \file cgnsWrite.cpp
 * \brief Main subroutines for cgns file write-out.
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

#ifdef USES_CGNS
#include "cgnsWrite.hpp"
#include "globalMesh.hpp"

#ifndef checkCGNSError
#define checkCGNSError(val)                 \
    if ((val) != CG_OK) {                   \
        const auto errMsg = cg_get_error(); \
        LFatal("CGNS error: %s", errMsg);   \
    }
#endif // !checkCUDAError

namespace OpenHurricane {
    createClassNameStr(cgnsWrite, "cgns");
    registerObjFty(solutionWrite, cgnsWrite, controller);
} // namespace OpenHurricane

OpenHurricane::cgnsWrite::cgnsWrite(const flowModel &flows, const iteration &iter,
                                    const runtimeMesh &mesh, const controller &cont)
    : solutionWrite(flows, iter, mesh, cont), fileType_(cgnsIO::IS_CGNS_HDF5_FILE),
      firstCallCGBC_(true), BCParentElem_(), BCParentElemPos_() {}

void OpenHurricane::cgnsWrite::writeToFile() const {
    const auto outN = outFile(".cgns");
    PLInfo("    Writting result to file: %s\n", outN.c_str());

    const auto &fzl = mesh().faceZones();
    integer offset = 0;
    if (HurMPI::parRun()) {
        offset = mesh().allMeshCellNumber();
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                continue;
            } else {
                integer csize = fzl[fzi].size();
                HurMPI::reduce(csize, MPI_SUM);
                if (HurMPI::master()) {
                    offset += csize;
                }
            }
        }
    } else {
        offset = mesh().allMeshCellNumber();
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                continue;
            } else {
                offset += fzl[fzi].size();
            }
        }
    }

    cgnsIO cff(outN, cgnsIO::WRITE);
    cff.setCGNSFileType(fileType_);
    int iB = 0;
    int iZ = 0;
    if (HurMPI::master()) {
        cff.open();
        iB = cff.writeBase("Base");
        iZ = cff.writeZone(iB, "Zone", mesh().allMeshNodeNumber(), offset, 0);
        cff.writeUnits(iB);
        if (iter().isSteadyFlow()) {
            checkCGNSError(
                cg_simulation_type_write(cff.cgnsFN(), iB, SimulationType_t::NonTimeAccurate));

            checkCGNSError(cg_biter_write(cff.cgnsFN(), iB, "BaseIterativeData", 1));
            checkCGNSError(cg_goto(cff.cgnsFN(), iB, "BaseIterativeData", 0, "end"));
            cgsize_t ii = 1;
            real time = 0;
            checkCGNSError(
                cg_array_write("TimeValues", cff.getDataTypeFromTheProgram(), 1, &ii, &time));
            checkCGNSError(cg_ziter_write(cff.cgnsFN(), iB, iZ, "ZoneIterativeData"));
            checkCGNSError(cg_goto(cff.cgnsFN(), iB, "Zone_t", iZ, "ZoneIterativeData", 0, "end"));
            char fsn[32] = "FlowSolution\0";
            cgsize_t ij[2] = {32, 1};
            checkCGNSError(
                cg_array_write("FlowSolutionPointers", DataType_t::Character, 2, ij, fsn));

        } else {
            checkCGNSError(
                cg_simulation_type_write(cff.cgnsFN(), iB, SimulationType_t::TimeAccurate));

            checkCGNSError(cg_biter_write(cff.cgnsFN(), iB, "BaseIterativeData", 1));
            checkCGNSError(cg_goto(cff.cgnsFN(), iB, "BaseIterativeData", 0, "end"));
            cgsize_t ii = 1;
            real time = iter().pTStep().totalTime();
            checkCGNSError(
                cg_array_write("TimeValues", cff.getDataTypeFromTheProgram(), 1, &ii, &time));
            checkCGNSError(cg_ziter_write(cff.cgnsFN(), iB, iZ, "ZoneIterativeData"));
            checkCGNSError(cg_goto(cff.cgnsFN(), iB, "Zone_t", iZ, "ZoneIterativeData", 0, "end"));
            char fsn[32] = "FlowSolution\0";

            checkCGNSError(
                cg_array_write("FlowSolutionPointers", DataType_t::Character, 2, &ii, fsn));
        }
    }
    writeCells(cff, iB, iZ);
    writeBCFaces(cff, iB, iZ);
    writePoints(cff, iB, iZ);
    writeSolutions(cff, iB, iZ);
    cff.close();
}
void OpenHurricane::cgnsWrite::writePoints(cgnsIO &cff, int B, int Z) const {
    if (HurMPI::parRun()) {
        List<pointField> tpp(mesh().pointZones().size());

        for (integer pzid = 0; pzid < mesh().pointZones().size(); ++pzid) {
            mesh().globalMeshInfo().globalPointIndeces().getNodeList(pzid, tpp[pzid]);
        }
        if (HurMPI::master()) {
            Array<real> components(mesh().allMeshNodeNumber());
            for (int j = 0; j < feature<point>::nElements_; ++j) {
                if (j == 0) {
                    integer nn = 0;
                    for (integer i = 0; i < tpp.size(); ++i) {
                        for (integer n = 0; n < tpp[i].size(); ++n) {
                            components[nn++] = tpp[i][n].x();
                        }
                    }
                    int xc;

                    checkCGNSError(cg_coord_write(cff.cgnsFN(), B, Z,
                                                  cff.getDataTypeFromTheProgram(), "CoordinateX",
                                                  components.data(), &xc));
                } else if (j == 1) {
                    integer nn = 0;
                    for (integer i = 0; i < tpp.size(); ++i) {
                        for (integer n = 0; n < tpp[i].size(); ++n) {
                            components[nn++] = tpp[i][n].y();
                        }
                    }
                    int yc;

                    checkCGNSError(cg_coord_write(cff.cgnsFN(), B, Z,
                                                  cff.getDataTypeFromTheProgram(), "CoordinateY",
                                                  components.data(), &yc));
                } else if (j == 2) {
                    integer nn = 0;
                    for (integer i = 0; i < tpp.size(); ++i) {
                        for (integer n = 0; n < tpp[i].size(); ++n) {
                            components[nn++] = tpp[i][n].z();
                        }
                    }
                    int zc;

                    checkCGNSError(cg_coord_write(cff.cgnsFN(), B, Z,
                                                  cff.getDataTypeFromTheProgram(), "CoordinateZ",
                                                  components.data(), &zc));
                }
            }
        }
    } else {
        const auto &ps = mesh().points();

        if (feature<point>::nElements_ == 3) {
            Array<real> componentsx = ps.component(0);
            int xc;

            checkCGNSError(cg_coord_write(cff.cgnsFN(), B, Z, cff.getDataTypeFromTheProgram(),
                                          "CoordinateX", componentsx.data(), &xc));
            Array<real> componentsy = ps.component(1);
            int yc;
            checkCGNSError(cg_coord_write(cff.cgnsFN(), B, Z, cff.getDataTypeFromTheProgram(),
                                          "CoordinateY", componentsy.data(), &yc));

            Array<real> componentsz = ps.component(2);
            int zc;
            checkCGNSError(cg_coord_write(cff.cgnsFN(), B, Z, cff.getDataTypeFromTheProgram(),
                                          "CoordinateZ", componentsz.data(), &zc));

        } else if (feature<point>::nElements_ == 2) {
            Array<real> componentsx = ps.component(0);
            int xc;
            checkCGNSError(cg_coord_write(cff.cgnsFN(), B, Z, cff.getDataTypeFromTheProgram(),
                                          "CoordinateX", componentsx.data(), &xc));
            Array<real> componentsy = ps.component(1);
            int yc;
            checkCGNSError(cg_coord_write(cff.cgnsFN(), B, Z, cff.getDataTypeFromTheProgram(),
                                          "CoordinateY", componentsy.data(), &yc));

        } else {
            LFatal("Unsupport point type with %d elements", feature<point>::nElements_);
        }
    }
}
void OpenHurricane::cgnsWrite::writeCells(cgnsIO &cff, int B, int Z) const {
    if (HurMPI::parRun()) {
        const auto &czl = mesh().cellZones();
        const auto &cells = mesh().cells();
        integer offset = 0;
        for (integer czi = 0; czi < czl.size(); ++czi) {
            const auto &cz = czl[czi];
            const auto shapeType = cff.convertElementType(cz.shapeType());
            integer csize = cz.size();
            HurMPI::reduce(csize, MPI_SUM);
            cgsize_t start = 0;
            cgsize_t end = 0;
            if (HurMPI::master()) {
                start = static_cast<cgsize_t>(offset + 1);
                end = static_cast<cgsize_t>(offset + csize);
                offset += csize;
            }
            int npe;
            int ntotal = 0;
            ElementType_t shapeType1;
            switch (shapeType) {
            case TETRA_4:
            case PYRA_5:
            case PENTA_6:
            case HEXA_8:
            case TRI_3:
            case QUAD_4:
                npe = cff.numNodes(shapeType);
                ntotal = cz.size() * npe;
                break;
            case MIXED:
                for (integer icell = cz.firstIndex(); icell <= cz.lastIndex(); ++icell) {
                    shapeType1 = cff.convertElementType(cells[icell].shapeType());
                    npe = cff.numNodes(shapeType1);
                    ntotal += (npe + 1);
                }
                break;
            case NGON_n:
            default:
                LFatal("Unsupported element type: %s", ElementTypeName[shapeType]);
            }
            integerList nSizeL;
            integerList displs;
            integer allSize = 0;
            nSizeL.resize(HurMPI::getProcSize(), Zero);
            nSizeL[HurMPI::getProcRank()] = ntotal;
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
            integerList elements(ntotal);
            ntotal = 0;
            integer icount = 0;
            switch (shapeType) {
            case TETRA_4:
            case PYRA_5:
            case PENTA_6:
            case HEXA_8:
            case TRI_3:
            case QUAD_4:
                npe = cff.numNodes(shapeType);
                ntotal = cz.size() * npe;
                for (integer icell = cz.firstIndex(); icell <= cz.lastIndex(); ++icell) {
                    for (integer i = 0; i < cells[icell].nodeSize(); ++i) {
                        elements[icount * npe + i] =
                            mesh().globalMeshInfo().globalPointIndeces().toGlobalIndex(
                                cells[icell].nodei(i)) +
                            1;
                    }
                    icount++;
                }
                break;
            case MIXED:
                for (integer icell = cz.firstIndex(); icell <= cz.lastIndex(); ++icell) {
                    shapeType1 = cff.convertElementType(cells[icell].shapeType());
                    npe = cff.numNodes(shapeType1);
                    elements[ntotal] = shapeType1;
                    for (integer i = 0; i < cells[icell].nodeSize(); ++i) {
                        elements[ntotal + i + 1] =
                            mesh().globalMeshInfo().globalPointIndeces().toGlobalIndex(
                                cells[icell].nodei(i)) +
                            1;
                    }
                    ntotal += (npe + 1);
                }
                break;
            case NGON_n:
            default:
                break;
            }
            integerArray allNNL;
            if (HurMPI::master()) {
                allNNL.resize(allSize);
            }
            HurMPI::Request request;
            HurMPI::igatherv(elements.data(), ntotal, feature<integer>::MPIType, allNNL.data(),
                             nSizeL.data(), displs.data(), feature<integer>::MPIType,
                             HurMPI::masterNo(), HurMPI::getComm(), &request);
            HurMPI::wait(&request, MPI_STATUSES_IGNORE);

            if (HurMPI::master()) {
                cgsize_t *allelements = nullptr;
                allelements = new cgsize_t[allSize];
                for (integer i = 0; i < allSize; ++i) {
                    allelements[i] = allNNL[i];
                }

                int SS;
                checkCGNSError(cg_section_write(cff.cgnsFN(), B, Z, cz.name().c_str(), shapeType,
                                                start, end, 0, allelements, &SS));

                HurDeleteDynArray(allelements);
            }
        }
    } else {
        const auto &czl = mesh().cellZones();
        const auto &cells = mesh().cells();
        integer nCell = 1;
        for (integer czi = 0; czi < czl.size(); ++czi) {
            const auto &cz = czl[czi];
            const auto shapeType = cff.convertElementType(cz.shapeType());

            auto start = cgsize_t(nCell + cz.firstIndex());
            auto end = cgsize_t(nCell + cz.lastIndex());
            nCell += cz.size();

            cgsize_t *elements = nullptr;
            int npe;
            int ntotal = 0;
            ElementType_t shapeType1;
            switch (shapeType) {
            case TETRA_4:
            case PYRA_5:
            case PENTA_6:
            case HEXA_8:
            case TRI_3:
            case QUAD_4:
                npe = cff.numNodes(shapeType);
                ntotal = cz.size() * npe;
                break;
            case MIXED:
                for (integer icell = cz.firstIndex(); icell <= cz.lastIndex(); ++icell) {
                    shapeType1 = cff.convertElementType(cells[icell].shapeType());
                    npe = cff.numNodes(shapeType1);
                    ntotal += (npe + 1);
                }
                break;
            case NGON_n:
            default:
                LFatal("Unsupported element type: %s", ElementTypeName[shapeType]);
            }
            elements = new cgsize_t[ntotal];
            ntotal = 0;
            integer icount = 0;
            switch (shapeType) {
            case TETRA_4:
            case PYRA_5:
            case PENTA_6:
            case HEXA_8:
            case TRI_3:
            case QUAD_4:
                npe = cff.numNodes(shapeType);
                for (integer icell = cz.firstIndex(); icell <= cz.lastIndex(); ++icell) {
                    for (integer i = 0; i < cells[icell].nodeSize(); ++i) {
                        elements[icount * npe + i] =
                            static_cast<cgsize_t>(cells[icell].nodei(i)) + 1;
                    }
                    icount++;
                }
                break;
            case MIXED:
                for (integer icell = cz.firstIndex(); icell <= cz.lastIndex(); ++icell) {
                    shapeType1 = cff.convertElementType(cells[icell].shapeType());
                    npe = cff.numNodes(shapeType1);
                    elements[ntotal] = shapeType1;
                    for (integer i = 0; i < cells[icell].nodeSize(); ++i) {
                        elements[ntotal + i + 1] = static_cast<cgsize_t>(cells[icell].nodei(i) + 1);
                    }
                    ntotal += (npe + 1);
                }
                break;
            case NGON_n:
            default:
                break;
            }
            int SS;
            checkCGNSError(cg_section_write(cff.cgnsFN(), B, Z, cz.name().c_str(), shapeType, start,
                                            end, 0, elements, &SS));
            HurDeleteDynArray(elements);
        }
    }
}
void OpenHurricane::cgnsWrite::writeBCFaces(cgnsIO &cff, int B, int Z) const {
    if (firstCallCGBC_) {
        firstCallCGBC_ = false;
        const auto &fzl = mesh().faceZones();
        const auto &faces = mesh().faces();
        const auto &cells = mesh().cells();
        integer countBC = 0;
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                continue;
            } else {
                countBC++;
            }
        }
        if (HurMPI::master()) {
            BCParentElem_.resize(countBC);
            BCParentElemPos_.resize(countBC);
        }

        if (HurMPI::parRun()) {
            countBC = 0;
            for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
                if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                    continue;
                } else {
                    integer csize = fzl[fzi].size();
                    int ntotal = csize * 2;

                    integerList nSizeL;
                    integerList displs;
                    integer allSize = 0;
                    nSizeL.resize(HurMPI::getProcSize(), Zero);
                    nSizeL[HurMPI::getProcRank()] = ntotal;
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

                    integerList elements(ntotal);
                    integer icount = 0;
                    for (integer iface = fzl[fzi].firstIndex(); iface <= fzl[fzi].lastIndex();
                         ++iface) {
                        elements[icount * 2 + 0] =
                            mesh().globalMeshInfo().globalCellIndeces().toGlobalIndex(
                                faces[iface].leftCell());
                        const auto cl = faces[iface].leftCell();

                        elements[icount * 2 + 1] = cells[cl].cgnsFaceConoId(iface);
                        icount++;
                    }
                    integerArray allNNL;
                    if (HurMPI::master()) {
                        allNNL.resize(allSize);
                        BCParentElem_[countBC].resize(allSize / 2);
                        BCParentElemPos_[countBC].resize(allSize / 2);
                    }
                    HurMPI::Request request;
                    HurMPI::igatherv(elements.data(), ntotal, feature<integer>::MPIType,
                                     allNNL.data(), nSizeL.data(), displs.data(),
                                     feature<integer>::MPIType, HurMPI::masterNo(),
                                     HurMPI::getComm(), &request);
                    HurMPI::wait(&request, MPI_STATUSES_IGNORE);
                    if (HurMPI::master()) {

                        for (integer i = 0; i < allSize / 2; ++i) {
                            BCParentElem_[countBC][i] = allNNL[i * 2 + 0];
                            BCParentElemPos_[countBC][i] = allNNL[i * 2 + 1];
                        }
                    }
                    countBC++;
                }
            }
        } else {
            countBC = 0;
            for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
                if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                    continue;
                } else {
                    integer csize = fzl[fzi].size();
                    BCParentElem_[countBC].resize(csize);
                    BCParentElemPos_[countBC].resize(csize);
                    integer icount = 0;
                    for (integer iface = fzl[fzi].firstIndex(); iface <= fzl[fzi].lastIndex();
                         ++iface) {
                        BCParentElem_[countBC][icount] = faces[iface].leftCell();
                        const auto cl = faces[iface].leftCell();
                        BCParentElemPos_[countBC][icount] = cells[cl].cgnsFaceConoId(iface);
                        icount++;
                    }
                    countBC++;
                }
            }
        }
    }
    if (HurMPI::parRun()) {
        const auto &fzl = mesh().faceZones();
        const auto &faces = mesh().faces();
        integer offset = mesh().allMeshCellNumber();
        integer countBC = 0;
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                continue;
            } else {
                integer csize = fzl[fzi].size();
                HurMPI::reduce(csize, MPI_SUM);
                cgsize_t start = 0;
                cgsize_t end = 0;
                cgsize_t pnts[2];
                if (HurMPI::master()) {
                    pnts[0] = static_cast<cgsize_t>(offset + 1);
                    pnts[1] = static_cast<cgsize_t>(offset + csize);
                    offset += csize;
                    int cc;
                    checkCGNSError(cg_boco_write(cff.cgnsFN(), B, Z, fzl[fzi].name().c_str(),
                                                 BCType_t::BCTypeUserDefined,
                                                 PointSetType_t::PointRange, 2, pnts, &cc));

                    checkCGNSError(cg_boco_gridlocation_write(cff.cgnsFN(), B, Z, cc,
                                                              GridLocation_t::FaceCenter));
                }
                const auto shapeType = cff.convertFaceElementType(fzl[fzi].faceType());
                int npe;
                int ntotal = 0;
                ElementType_t shapeType1;
                switch (shapeType) {
                case TRI_3:
                case QUAD_4:
                case HEXA_8:
                case BAR_2:
                    npe = cff.numNodes(shapeType);
                    ntotal = fzl[fzi].size() * npe;
                    break;
                case MIXED:
                    for (integer iface = fzl[fzi].firstIndex(); iface <= fzl[fzi].lastIndex();
                         ++iface) {
                        if (faces[iface].isTriangular()) {
                            shapeType1 = ElementType_t::TRI_3;
                        } else if (faces[iface].isQuadrilateral()) {
                            shapeType1 = ElementType_t::QUAD_4;
                        } else if (faces[iface].isPolygonal()) {
                            shapeType1 = ElementType_t::NFACE_n;
                        } else if (faces[iface].isLinear()) {
                            shapeType1 = ElementType_t::BAR_2;
                        } else {
                            LFatal("Unsupported element type: %s", ElementTypeName[shapeType1]);
                        }
                        npe = cff.numNodes(shapeType1);
                        ntotal += (npe + 1);
                    }
                    break;
                case NFACE_n:
                default:
                    LFatal("Unsupported element type: %s", ElementTypeName[shapeType]);
                }
                integerList nSizeL;
                integerList displs;
                integer allSize = 0;
                nSizeL.resize(HurMPI::getProcSize(), Zero);
                nSizeL[HurMPI::getProcRank()] = ntotal;
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
                integerList elements(ntotal);
                ntotal = 0;
                integer icount = 0;
                switch (shapeType) {
                case TRI_3:
                case QUAD_4:
                case HEXA_8:
                case BAR_2:
                    npe = cff.numNodes(shapeType);
                    ntotal = fzl[fzi].size() * npe;
                    for (integer iface = fzl[fzi].firstIndex(); iface <= fzl[fzi].lastIndex();
                         ++iface) {
                        for (integer i = 0; i < faces[iface].size(); ++i) {
                            elements[icount * npe + i] =
                                mesh().globalMeshInfo().globalPointIndeces().toGlobalIndex(
                                    faces[iface][i]) +
                                1;
                        }
                        icount++;
                    }
                    break;
                case MIXED:
                    for (integer iface = fzl[fzi].firstIndex(); iface <= fzl[fzi].lastIndex();
                         ++iface) {
                        if (faces[iface].isTriangular()) {
                            shapeType1 = ElementType_t::TRI_3;
                        } else if (faces[iface].isQuadrilateral()) {
                            shapeType1 = ElementType_t::QUAD_4;
                        } else if (faces[iface].isPolygonal()) {
                            shapeType1 = ElementType_t::NFACE_n;
                        } else if (faces[iface].isLinear()) {
                            shapeType1 = ElementType_t::BAR_2;
                        } else {
                            LFatal("Unsupported element type: %s", ElementTypeName[shapeType1]);
                        }
                        npe = cff.numNodes(shapeType1);
                        elements[ntotal] = shapeType1;
                        for (integer i = 0; i < faces[iface].size(); ++i) {
                            elements[ntotal + i + 1] =
                                mesh().globalMeshInfo().globalPointIndeces().toGlobalIndex(
                                    faces[iface][i] + 1);
                        }
                        ntotal += (npe + 1);
                    }
                    break;
                case NFACE_n:
                default:
                    break;
                }
                integerArray allNNL;
                if (HurMPI::master()) {
                    allNNL.resize(allSize);
                }
                HurMPI::Request request;
                HurMPI::igatherv(elements.data(), ntotal, feature<integer>::MPIType, allNNL.data(),
                                 nSizeL.data(), displs.data(), feature<integer>::MPIType,
                                 HurMPI::masterNo(), HurMPI::getComm(), &request);
                HurMPI::wait(&request, MPI_STATUSES_IGNORE);
                if (HurMPI::master()) {
                    cgsize_t *allelements = nullptr;
                    allelements = new cgsize_t[allSize];
                    for (integer i = 0; i < allSize; ++i) {
                        allelements[i] = allNNL[i];
                    }
                    int SS;
                    checkCGNSError(cg_section_write(cff.cgnsFN(), B, Z, fzl[fzi].name().c_str(),
                                                    shapeType, pnts[0], pnts[1], 0, allelements,
                                                    &SS));

                    HurDeleteDynArray(allelements);

                    List<cgsize_t> parentData(BCParentElem_[countBC].size() * 4, cgsize_t(0));
                    for (integer i = 0; i < BCParentElem_[countBC].size(); ++i) {
                        parentData[i] = static_cast<cgsize_t>(BCParentElem_[countBC][i] + 1);
                        parentData[BCParentElem_[countBC].size() * 2 + i] =
                            static_cast<cgsize_t>(BCParentElemPos_[countBC][i] + 1);
                    }
                    checkCGNSError(cg_parent_data_write(cff.cgnsFN(), B, Z, SS, parentData.data()));
                }
                countBC++;
            }
        }
    } else {
        const auto &fzl = mesh().faceZones();
        const auto &faces = mesh().faces();
        integer offset = mesh().allMeshCellNumber();
        integer countBC = 0;
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                continue;
            } else {
                cgsize_t pnts[2];
                pnts[0] = static_cast<cgsize_t>(offset + 1);
                pnts[1] = static_cast<cgsize_t>(fzl[fzi].size() + offset);
                offset += fzl[fzi].size();
                int cc;
                checkCGNSError(cg_boco_write(cff.cgnsFN(), B, Z, fzl[fzi].name().c_str(),
                                             BCType_t::BCTypeUserDefined,
                                             PointSetType_t::PointRange, 2, pnts, &cc));
                checkCGNSError(
                    cg_boco_gridlocation_write(cff.cgnsFN(), B, Z, cc, GridLocation_t::FaceCenter));

                const auto shapeType = cff.convertFaceElementType(fzl[fzi].faceType());
                cgsize_t *elements = nullptr;
                int npe;
                int ntotal = 0;
                ElementType_t shapeType1 = ElementTypeNull;
                switch (shapeType) {
                case TRI_3:
                case QUAD_4:
                case HEXA_8:
                case BAR_2:
                    npe = cff.numNodes(shapeType);
                    ntotal = fzl[fzi].size() * npe;
                    break;
                case MIXED:
                    for (integer iface = fzl[fzi].firstIndex(); iface <= fzl[fzi].lastIndex();
                         ++iface) {
                        if (faces[iface].isTriangular()) {
                            shapeType1 = ElementType_t::TRI_3;
                        } else if (faces[iface].isQuadrilateral()) {
                            shapeType1 = ElementType_t::QUAD_4;
                        } else if (faces[iface].isPolygonal()) {
                            shapeType1 = ElementType_t::NFACE_n;
                        } else if (faces[iface].isLinear()) {
                            shapeType1 = ElementType_t::BAR_2;
                        } else {
                            LFatal("Unsupported element type: %s", ElementTypeName[shapeType1]);
                        }
                        npe = cff.numNodes(shapeType1);
                        ntotal += (npe + 1);
                    }
                    break;
                case NFACE_n:
                default:
                    LFatal("Unsupported element type: %s", ElementTypeName[shapeType]);
                }
                elements = new cgsize_t[ntotal];
                ntotal = 0;
                integer icount = 0;
                switch (shapeType) {
                case TRI_3:
                case QUAD_4:
                case HEXA_8:
                case BAR_2:
                    npe = cff.numNodes(shapeType);
                    for (integer iface = fzl[fzi].firstIndex(); iface <= fzl[fzi].lastIndex();
                         ++iface) {
                        for (integer i = 0; i < faces[iface].size(); ++i) {
                            elements[icount * npe + i] = static_cast<cgsize_t>(faces[iface][i]) + 1;
                        }
                        icount++;
                    }
                    break;
                case MIXED:
                    for (integer iface = fzl[fzi].firstIndex(); iface <= fzl[fzi].lastIndex();
                         ++iface) {
                        if (faces[iface].isTriangular()) {
                            shapeType1 = ElementType_t::TRI_3;
                        } else if (faces[iface].isQuadrilateral()) {
                            shapeType1 = ElementType_t::QUAD_4;
                        } else if (faces[iface].isPolygonal()) {
                            shapeType1 = ElementType_t::NFACE_n;
                        } else if (faces[iface].isLinear()) {
                            shapeType1 = ElementType_t::BAR_2;
                        } else {
                            LFatal("Unsupported element type: %s", ElementTypeName[shapeType1]);
                        }
                        npe = cff.numNodes(shapeType1);
                        elements[ntotal] = shapeType1;
                        for (integer i = 0; i < faces[iface].size(); ++i) {
                            elements[ntotal + i + 1] = static_cast<cgsize_t>(faces[iface][i] + 1);
                        }
                        ntotal += (npe + 1);
                    }
                    break;
                case NFACE_n:
                default:
                    break;
                }
                int SS;
                checkCGNSError(cg_section_write(cff.cgnsFN(), B, Z, fzl[fzi].name().c_str(),
                                                shapeType, pnts[0], pnts[1], 0, elements, &SS));

                HurDeleteDynArray(elements);
                List<cgsize_t> parentData(BCParentElem_[countBC].size() * 4, cgsize_t(0));
                for (integer i = 0; i < BCParentElem_[countBC].size(); ++i) {
                    parentData[i] = static_cast<cgsize_t>(BCParentElem_[countBC][i] + 1);
                    parentData[BCParentElem_[countBC].size() * 2 + i] =
                        static_cast<cgsize_t>(BCParentElemPos_[countBC][i] + 1);
                }
                checkCGNSError(cg_parent_data_write(cff.cgnsFN(), B, Z, SS, parentData.data()));
                countBC++;
            }
        }
    }
}
void OpenHurricane::cgnsWrite::writeSolutions(cgnsIO &cff, int B, int Z) const {
    int SS = 0;
    if (HurMPI::master()) {
        checkCGNSError(
            cg_sol_write(cff.cgnsFN(), B, Z, "FlowSolution", GridLocation_t ::CellCenter, &SS));
    }

    auto orm = mesh().outResultMap(this->outVarNameList());
    auto outf = [this, &cff, B, Z, SS](object *obj) {
        const auto &on = obj->outputVarNameL();
        for (int d = 0; d < obj->nElements(); ++d) {
            auto f = obj->realComponent(d);
            auto FF = writeArray(cff, B, Z, SS, f, on[d].c_str());
        }
    };
    for (integer i = 0; i < orm.size(); ++i) {
        outf(orm[i]);
    }
}

int OpenHurricane::cgnsWrite::writeArray(cgnsIO &cff, int B, int Z, int S, realArray &tmpv,
                                         const char *nstr) const {
    if (HurMPI::parRun()) {
        const auto &fzl = mesh().faceZones();
        const auto &faces = mesh().faces();
        integer offset = mesh().allMeshCellNumber();
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                continue;
            } else {
                integer csize = fzl[fzi].size();
                HurMPI::reduce(csize, MPI_SUM);
                if (HurMPI::master()) {
                    offset += csize;
                }
                const auto &fz = fzl[fzi];
                for (integer iface = fz.firstIndex(); iface <= fz.lastIndex(); ++iface) {
                    const auto cl = faces[iface].leftCell();
                    const auto cr = faces[iface].rightCell();
                    tmpv[cr] = 0.5 * (tmpv[cl] + tmpv[cr]);
                }
            }
        }

        integerList nSizeL;
        integerList displs;
        integer allSize = 0;
        nSizeL.resize(HurMPI::getProcSize(), Zero);
        nSizeL[HurMPI::getProcRank()] = mesh().nCells();
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
        realArray allTmpv;

        if (HurMPI::master()) {
            allTmpv.resize(offset);
        }
        HurMPI::Request request;
        HurMPI::igatherv(tmpv.data(), mesh().nCells(), feature<real>::MPIType, allTmpv.data(),
                         nSizeL.data(), displs.data(), feature<real>::MPIType, HurMPI::masterNo(),
                         HurMPI::getComm(), &request);
        HurMPI::wait(&request, MPI_STATUSES_IGNORE);
        offset = mesh().allMeshCellNumber();
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                continue;
            } else {
                integer csize = fzl[fzi].size();
                const auto &fz = fzl[fzi];
                realArray tmpvBc(csize);
                integer count = 0;
                for (integer iface = fz.firstIndex(); iface <= fz.lastIndex(); ++iface) {
                    const auto cr = faces[iface].rightCell();
                    tmpvBc[count++] = tmpv[cr];
                }

                integerList nSizeL0;
                integerList displs0;
                integer allSize0 = 0;
                nSizeL0.resize(HurMPI::getProcSize(), Zero);
                nSizeL0[HurMPI::getProcRank()] = csize;
                if (HurMPI::master()) {
                    displs0.resize(HurMPI::getProcSize(), Zero);
                }
                HurMPI::gatherList(nSizeL0, HurMPI::masterNo(), HurMPI::getComm());
                if (HurMPI::master()) {
                    for (integer ip = 1; ip < HurMPI::getProcSize(); ++ip) {
                        displs0[ip] = displs0[ip - 1] + nSizeL0[ip - 1];
                    }
                    for (integer ip = 0; ip < HurMPI::getProcSize(); ++ip) {
                        allSize0 += nSizeL0[ip];
                    }
                }
                realArray allTmpvbc;

                if (HurMPI::master()) {
                    allTmpvbc.resize(allSize0);
                }
                HurMPI::Request request0;
                HurMPI::igatherv(tmpvBc.data(), csize, feature<real>::MPIType, allTmpvbc.data(),
                                 nSizeL0.data(), displs0.data(), feature<real>::MPIType,
                                 HurMPI::masterNo(), HurMPI::getComm(), &request0);
                HurMPI::wait(&request0, MPI_STATUSES_IGNORE);

                HurMPI::reduce(csize, MPI_SUM);
                for (integer i = 0; i < allSize0; ++i) {
                    allTmpv[offset + i] = allTmpvbc[i];
                }
                if (HurMPI::master()) {
                    offset += csize;
                }
            }
        }

        int FF = 0;
        if (HurMPI::master()) {
            checkCGNSError(cg_field_write(cff.cgnsFN(), B, Z, S, cff.getDataTypeFromTheProgram(),
                                          nstr, allTmpv.data(), &FF));
        }
        return FF;

    } else {
        const auto &fzl = mesh().faceZones();
        const auto &faces = mesh().faces();
        integer countBC = 0;
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                continue;
            } else {
                const auto &fz = fzl[fzi];
                countBC += fz.size();
            }
        }

        realArray tmpvbc(countBC);
        countBC = 0;
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                continue;
            } else {
                const auto &fz = fzl[fzi];
                for (integer iface = fz.firstIndex(); iface <= fz.lastIndex(); ++iface) {
                    const auto cl = faces[iface].leftCell();
                    const auto cr = faces[iface].rightCell();
                    tmpvbc[countBC++] = 0.5 * (tmpv[cl] + tmpv[cr]);
                }
            }
        }
        countBC = 0;
        for (integer fzi = 0; fzi < fzl.size(); ++fzi) {
            if (fzl[fzi].isInterior() || fzl[fzi].isCutFace()) {
                continue;
            } else {
                const auto &fz = fzl[fzi];
                for (integer iface = fz.firstIndex(); iface <= fz.lastIndex(); ++iface) {

                    tmpv[countBC + mesh().nCells()] = tmpvbc[countBC];
                    countBC++;
                }
            }
        }
        int FF;
        checkCGNSError(cg_field_write(cff.cgnsFN(), B, Z, S, cff.getDataTypeFromTheProgram(), nstr,
                                      tmpv.data(), &FF));
        return FF;
    }
    return 0;
}
#endif // USES_CGNS