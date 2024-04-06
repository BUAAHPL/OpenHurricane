/*!
 * \file cgnsIOWrite.cpp
 * \brief Main subroutines of the <i>cgnsIOWrite.hpp</i> file.
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

int OpenHurricane::cgnsIO::writeBase(const string &baseName, const integer cellDim,
                                 const integer phyDim) {
    if (closed()) {
        errorAbortStr(
            ("Attempt to write base information: " + baseName + " to closed file: " + filename_));
    }
    if (isRead()) {
        errorAbortStr(("Attempt to write base information: " + baseName +
                       " to read only file: " + filename_));
    }
    int baseIndex = -1;
    const auto result = cg_base_write(cgnsFN_, baseName.c_str(), cellDim, phyDim, &baseIndex);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return baseIndex;
}

int OpenHurricane::cgnsIO::writeZone(const int B, const string &zoneName, cgsize_t *size,
                                 ZoneType_t zoneType) {
    if (closed()) {
        errorAbortStr(
            ("Attempt to write zone information: " + zoneName + " to closed file: " + filename_));
    }
    if (isRead()) {
        errorAbortStr(("Attempt to write zone information: " + zoneName +
                       " to read only file: " + filename_));
    }
    int zoneIndex = -1;
    const auto result = cg_zone_write(cgnsFN_, B, zoneName.c_str(), size, zoneType, &zoneIndex);
    if (result != CG_OK) {
        const auto errMsg = cg_get_error();
        LFatal(errMsg);
    }
    return zoneIndex;
}

int OpenHurricane::cgnsIO::writeZone(const int B, const string &zoneName, const integer &NVertex,
                                 const integer &NCell3D, const integer &NBoundVertex) {
    cgsize_t *size = new cgsize_t[3];
    size[0] = NVertex;
    size[1] = NCell3D;
    size[2] = NBoundVertex;
    const int Z = writeZone(B, zoneName, size, Unstructured);
    delete[] size;
    return Z;
}

void OpenHurricane::cgnsIO::writeSimulationType(const int B, const SimulationType_t st) {
    if (closed()) {
        LFatal("Attempt to write simulation type to closed file: %s", filename_.c_str());
    }
    if (isRead()) {
        LFatal("Attempt to write simulation type to read only file: %s", filename_.c_str());
    }
    const auto result = cg_simulation_type_write(cgnsFN_, B, st);
    if (result != CG_OK) {
        LFatal(cg_get_error());
    }
}

void OpenHurricane::cgnsIO::writeBndConLocation(const int B, const int Z, const int BC,
                                            const CGNSGridLocation &gl) const {
    if (closed()) {
        LFatal("Attempt to write simulation type to closed file: %s", filename_.c_str());
    }
    if (isRead()) {
        LFatal("Attempt to write simulation type to read only file: %s", filename_.c_str());
    }
    const auto result = cg_boco_gridlocation_write(cgnsFN_, B, Z, BC, convertGridLocation(gl));
    if (result != CG_OK) {
        LFatal(cg_get_error());
    }
}

void OpenHurricane::cgnsIO::writeUnits(const int B) const {
    if (closed()) {
        LFatal("Attempt to write simulation type to closed file: %s", filename_.c_str());
    }
    if (isRead()) {
        LFatal("Attempt to write simulation type to read only file: %s", filename_.c_str());
    }
    auto result = cg_goto(cgnsFN_, B, "end");
    if (result != CG_OK) {
        LFatal(cg_get_error());
    }
    result = cg_units_write(MassUnits_t::Kilogram, LengthUnits_t::Meter, TimeUnits_t::Second,
                            TemperatureUnits_t::Kelvin, AngleUnits_t::Radian);
    if (result != CG_OK) {
        LFatal(cg_get_error());
    }
    result = cg_dataclass_write(DataClass_t::Dimensional);
    if (result != CG_OK) {
        LFatal(cg_get_error());
    }
}

#endif // USES_CGNS