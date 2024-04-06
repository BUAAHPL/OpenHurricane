/*!
 * \file cgnsIORead.cpp
 * \brief Main subroutines of the <i>cgnsIORead.hpp</i> file.
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
#include "cgnsIO.hpp"

#ifdef USES_CGNS

hur_nodiscard int OpenHurricane::cgnsIO::readNBases() const {
    if (closed()) {
        LFatal("Attempt to read the number of bases from closed file: %s", filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of bases from write only file: %s", filename_.c_str());
    }
    int baseIndex = -1;

    const auto result = cg_nbases(cgnsFN_, &baseIndex);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return baseIndex;
}

void OpenHurricane::cgnsIO::readBase(const int B, string &baseName, integer &cellDim,
                                 integer &phyDim) const {
    if (closed()) {
        LFatal("Attempt to read the base from closed file: %s", filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the base from write only file: %s", filename_.c_str());
    }
    char baseN[33];
    const auto result = cg_base_read(cgnsFN_, B, baseN, &cellDim, &phyDim);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    baseName = baseN;
}

void OpenHurricane::cgnsIO::readCellDim(const int B, integer &cellDim) const {
    if (closed()) {
        LFatal("Attempt to read the dimension of cells from closed file: %s", filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the dimension of cells from write only file: %s",
               filename_.c_str());
    }
    const auto result = cg_cell_dim(cgnsFN_, B, &cellDim);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
}

hur_nodiscard int OpenHurricane::cgnsIO::readNZones(const int B) const {
    if (closed()) {
        LFatal("Attempt to read the number of zones from closed file: %s", filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of zones from write only file: %s", filename_.c_str());
    }
    int nzones = -1;

    const auto result = cg_nzones(cgnsFN_, B, &nzones);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return nzones;
}

hur_nodiscard int OpenHurricane::cgnsIO::readIndexDim(const int B, const int Z) const {
    if (closed()) {
        LFatal("Attempt to read the index dimension of zone from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the index dimension of zone from write only file: %s",
               filename_.c_str());
    }
    int indexDim = -1;

    const auto result = cg_index_dim(cgnsFN_, B, Z, &indexDim);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return indexDim;
}

void OpenHurricane::cgnsIO::readZoneType(const int B, const int Z, ZoneType_t *zoneType) const {
    if (closed()) {
        LFatal("Attempt to read the zonetype from closed file: %s", filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the zonetype from write only file: %s", filename_.c_str());
    }
    const auto result = cg_zone_type(cgnsFN_, B, Z, zoneType);

    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
}

void OpenHurricane::cgnsIO::readZone(const int B, const int Z, string &zoneName, cgsize_t *size) const {
    if (closed()) {
        LFatal("Attempt to read the zone info from closed file: %s", filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the zone info from write only file: %s", filename_.c_str());
    }
    char zncPtr[33];
    const auto result = cg_zone_read(cgnsFN_, B, Z, zncPtr, size);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    zoneName = zncPtr;
}

void OpenHurricane::cgnsIO::readZone(const int B, const int Z, string &zoneName, integer &NVertex,
                                 integer &NCell3D, integer &NBoundVertex) const {
    if (!isUnstructuredZone(B, Z)) {
        LFatal("The current program does not support structured zone");
    }
    cgsize_t size[3];
    readZone(B, Z, zoneName, size);
    NVertex = (integer)size[0];
    NCell3D = (integer)size[1];
    NBoundVertex = (integer)size[2];
}

void OpenHurricane::cgnsIO::readSimulationType(const int B, SimulationType_t &st) const {
    if (closed()) {
        LFatal("Attempt to read the simulation type from closed file: %s", filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the simulation type from write only file: %s", filename_.c_str());
    }
    const auto result = cg_simulation_type_read(cgnsFN_, B, &st);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
}

hur_nodiscard int OpenHurricane::cgnsIO::readNCoords(const int B, const int Z) const {
    if (closed()) {
        LFatal("Attempt to read the number of coordinate arrays from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of coordinate arrays from write only file: %s",
               filename_.c_str());
    }
    int nc;
    const auto result = cg_ncoords(cgnsFN_, B, Z, &nc);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return nc;
}

void OpenHurricane::cgnsIO::readCoordInfo(const int B, const int Z, const int C, DataType_t &dt,
                                      string &coordName) const {
    if (closed()) {
        LFatal("Attempt to read the info of coordinate arrays from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the info of coordinate arrays from write only file: %s",
               filename_.c_str());
    }
    char cn[33];
    const auto result = cg_coord_info(cgnsFN_, B, Z, C, &dt, cn);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    coordName = cn;
}

void OpenHurricane::cgnsIO::readCoord(const int B, const int Z, const integer NVertex,
                                  const DataType_t dt, const string &coordName,
                                  realArray &coord) const {
    if (coord.size() != NVertex) {
        coord.resize(NVertex);
    }
    if (closed()) {
        LFatal("Attempt to read the info of coordinate arrays from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the info of coordinate arrays from write only file: %s",
               filename_.c_str());
    }
    cgsize_t range_min = 1, range_max = NVertex;
    const auto result =
        cg_coord_read(cgnsFN_, B, Z, coordName.c_str(), dt, &range_min, &range_max, coord.data());
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::cgnsIO::readCoord(const int B, const int Z,
                                                                const integer NVertex,
                                                                const DataType_t dt,
                                                                const string &coordName) const {
    realArray coord(NVertex);
    readCoord(B, Z, NVertex, dt, coordName, coord);
    return coord;
}

void OpenHurricane::cgnsIO::readCoord(const int B, const int Z, const integer NVertex,
                                  vectorArray &coord) const {
    integer nCA = readNCoords(B, Z);
    if (coord.size() != NVertex) {
        coord.resize(NVertex);
    }
    realArray coordComp(NVertex);
#ifdef HUR_DEBUG
    Pout("    Info: the number of coordinate arrays is %d", nCA);
#endif // HUR_DEBUG
    for (int kk = 0; kk < nCA; ++kk) {
        DataType_t dtt;
        string coordName;
        readCoordInfo(B, Z, kk + 1, dtt, coordName);
#ifdef HUR_DEBUG
        Pout << "    Info: the name of coordinate array " << kk + 1 << " is " << coordName
             << std::endl;
#endif // HUR_DEBUG
        readCoord(B, Z, NVertex, getDataTypeFromTheProgram(), coordName, coordComp);
        if (coordName == "CoordinateX") {
            coord.setComponent(0, coordComp);
        } else if (coordName == "CoordinateY") {
            coord.setComponent(1, coordComp);
        } else if (coordName == "CoordinateZ") {
            coord.setComponent(2, coordComp);
        } else {
            errorAbortStr(("Unsupported Coordinate Identifier: " + coordName));
        }
    }
}

hur_nodiscard int OpenHurricane::cgnsIO::readNSections(const int B, const int Z) const {
    if (closed()) {
        LFatal("Attempt to read the number of element sections from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of element sections from write only file: %s",
               filename_.c_str());
    }
    int nc;
    const auto result = cg_nsections(cgnsFN_, B, Z, &nc);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return nc;
}

void OpenHurricane::cgnsIO::readSection(const int B, const int Z, const int S, string &eleSectName,
                                    ElementType_t *type, integer &start, integer &end, int *nbndry,
                                    int *parentFlag) const {
    if (closed()) {
        LFatal("Attempt to read element sections from closed file: %s", filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read element sections from write only file: %s", filename_.c_str());
    }
    char eln[33];
    cgsize_t cgStart = 0, cgEnd = 0;
    const auto result =
        cg_section_read(cgnsFN_, B, Z, S, eln, type, &cgStart, &cgEnd, nbndry, parentFlag);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    eleSectName = eln;
    start = (integer)cgStart;
    end = (integer)cgEnd;
}

hur_nodiscard OpenHurricane::integer OpenHurricane::cgnsIO::readElementDataSize(const int B, const int Z,
                                                                        const int S) const {
    if (closed()) {
        LFatal("Attempt to read element sections from closed file: %s", filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read element sections from write only file: %s", filename_.c_str());
    }
    cgsize_t cgEs = 0;
    const auto result = cg_ElementDataSize(cgnsFN_, B, Z, S, &cgEs);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return integer(cgEs);
}

void OpenHurricane::cgnsIO::readElement(const int B, const int Z, const int S, cgsize_t *elementData,
                                    cgsize_t *parentData) const {
    if (closed()) {
        LFatal("Attempt to read element from closed file: %s", filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read element from write only file: %s", filename_.c_str());
    }
    cgsize_t cgEs = 0;
    const auto result = cg_elements_read(cgnsFN_, B, Z, S, elementData, parentData);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
}

hur_nodiscard int OpenHurricane::cgnsIO::readNBndCon(const int B, const int Z) const {
    if (closed()) {
        LFatal("Attempt to read the number of boundary conditions from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of boundary conditions from write only file: %s",
               filename_.c_str());
    }
    int nbc = 0;
    const auto result = cg_nbocos(cgnsFN_, B, Z, &nbc);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return nbc;
}

hur_nodiscard double OpenHurricane::cgnsIO::readBndConId(const int B, const int Z, const int BC) const {
    if (closed()) {
        LFatal("Attempt to read the number of boundary conditions from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of boundary conditions from write only file: %s",
               filename_.c_str());
    }
    double id;
    const auto result = cg_boco_id(cgnsFN_, B, Z, BC, &id);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return id;
}

hur_nodiscard OpenHurricane::cgnsIO::CGNSGridLocation
OpenHurricane::cgnsIO::readBndConLocation(const int B, const int Z, const int BC) const {
    if (closed()) {
        LFatal("Attempt to read the number of boundary conditions from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of boundary conditions from write only file: %s",
               filename_.c_str());
    }
    GridLocation_t gl;
    const auto result = cg_boco_gridlocation_read(cgnsFN_, B, Z, BC, &gl);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }

    return convertGridLocation(gl);
}

void OpenHurricane::cgnsIO::readBndConInfo(const int B, const int Z, const int BC,
                                       faceZoneList &fzl) const {
    if (closed()) {
        LFatal("Attempt to read the number of boundary conditions from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of boundary conditions from write only file: %s",
               filename_.c_str());
    }
    char boconame[33];
    BCType_t bocotype;
    PointSetType_t ptset_type;
    cgsize_t npnts;
    int NormalIndex;
    cgsize_t NormalListSize;
    DataType_t NormalDataType;
    int ndataset;
    const auto result = cg_boco_info(cgnsFN_, B, Z, BC, boconame, &bocotype, &ptset_type, &npnts,
                                     &NormalIndex, &NormalListSize, &NormalDataType, &ndataset);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }

    bool found = false;
    for (auto &e : fzl) {
        if (e.name() == boconame) {
            e.setBcType(checkBCType(bocotype));

            found = true;
            break;
        }
    }
    if (!found) {
        if (cg_goto(cgnsFN_, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BC, NULL) != CG_OK) {
            const auto errMsg = cg_get_error();
            LFatal(errMsg);
        }
        /*auto nnf = readNodeNFamilies();*/

        char familyame[33];
        if (cg_famname_read(familyame) != CG_OK) {
            const auto errMsg = cg_get_error();
            LFatal(errMsg);
        }
        for (auto &e : fzl) {
            if (e.name() == familyame) {
                e.setBcType(checkBCType(bocotype));
                found = true;
                break;
            }
        }
#ifdef HUR_DEBUG
        if (!found) {
            checkWarningStr(("Cann not find \"" + std::string(familyame) + "\" in face zones"));
        }
#endif // HUR_DEBUG
    }
}

hur_nodiscard int OpenHurricane::cgnsIO::readNConns(const int B, const int Z) const {
    if (closed()) {
        LFatal("Attempt to read the number of boundary conditions from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of boundary conditions from write only file: %s",
               filename_.c_str());
    }
    int nf;
    const auto result = cg_nconns(cgnsFN_, B, Z, &nf);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return nf;
}

void OpenHurricane::cgnsIO::readConnInfo(const int B, const int Z, const int I,
                                     string &connectName) const {
    if (closed()) {
        LFatal("Attempt to read the number of boundary conditions from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of boundary conditions from write only file: %s",
               filename_.c_str());
    }
    char connectname[33];
    GridLocation_t location;
    GridConnectivityType_t connect_type;
    PointSetType_t ptset_type;
    cgsize_t npnts;
    char donorname[33];
    ZoneType_t donor_zonetype;
    PointSetType_t donor_ptset_type;
    DataType_t donor_datatype;
    cgsize_t ndata_donor;
    const auto result =
        cg_conn_info(cgnsFN_, B, Z, I, connectname, &location, &connect_type, &ptset_type, &npnts,
                     donorname, &donor_zonetype, &donor_ptset_type, &donor_datatype, &ndata_donor);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    connectName = connectname;
}

hur_nodiscard int OpenHurricane::cgnsIO::readNFamilies(const int B) const {
    if (closed()) {
        LFatal("Attempt to read the number of boundary conditions from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of boundary conditions from write only file: %s",
               filename_.c_str());
    }
    int nf;
    const auto result = cg_nfamilies(cgnsFN_, B, &nf);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return nf;
}

hur_nodiscard int OpenHurricane::cgnsIO::readNodeNFamilies() const {
    if (closed()) {
        LFatal("Attempt to read the number of boundary conditions from closed file: %s",
               filename_.c_str());
    }
    if (isWrite()) {
        LFatal("Attempt to read the number of boundary conditions from write only file: %s",
               filename_.c_str());
    }
    int nf;
    const auto result = cg_node_nfamilies(&nf);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return nf;
}

#endif // USES_CGNS