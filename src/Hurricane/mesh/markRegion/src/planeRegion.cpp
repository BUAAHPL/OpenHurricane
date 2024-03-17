/*!
 * \file planeRegion.cpp
 * \brief Main subroutines for plane region of mesh.
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

#include "planeRegion.hpp"
#include "cellMesh.hpp"
#include "formulaParsing.hpp"
#include "constants.hpp"


namespace OpenHurricane {
    createClassNameStr(planeRegion, "plane");
}
namespace OpenHurricane {
    registerObjFty(markRegion, planeRegion, controller);
}

OpenHurricane::planeRegion::planeRegion(const controller &cont)
    : markRegion(cont), normal_(cont.findOrDefault<vector>("normal", vector(0.0))),
      planePoint_(cont.findOrDefault<vector>("point", vector(0.0))) {}

hur_nodiscard OpenHurricane::integerList
OpenHurricane::planeRegion::regionCellId(const runtimeMesh &mesh) const {
    integerList cid;
    if (!isOptionSet()) {
        PLWarning("Option unset for mark region: %d", id());
        return cid;
    }

    vector sc = normal_.normalized();
    vector pp = planePoint_;

    const auto &cC = mesh.cellCentre();
    cid.resize(mesh.nCells(), -1);
    integer count = 0;

    for (integer n = 0; n < mesh.nCells(); ++n) {
        real proj = (cC[n] - pp) * sc;

        if (isInside()) {
            if (proj > 0.0)
                continue;
        } else {
            if (proj <= 0.0)
                continue;
        }
        cid[count] = n;
        count++;
    }
    cid.resize(count);
    return cid;
}

bool OpenHurricane::planeRegion::patching(realGeometryArray<cellMesh> &cellQ, real &value) const {
    return patch(cellQ, value);
}

bool OpenHurricane::planeRegion::patching(vectorGeometryArray<cellMesh> &cellQ, vector &value) const {
    return patch(cellQ, value);
}

bool OpenHurricane::planeRegion::distributing(realGeometryArray<cellMesh> &cellQ,
                                              std::string &value) const {
    const auto &mesh = cellQ.mesh();
    std::string oriof = value;

    if (!isOptionSet()) {
        PLWarning("Option unset for mark region: %d", id());
        return false;
    }

    vector sc = normal_.normalized();
    vector pp = planePoint_;

    const auto &cC = mesh.cellCentre();
    stdStringList vf(3);

    for (integer n = 0; n < mesh.nCells(); ++n) {
        std::string of = oriof;
        real proj = (cC[n] - pp) * sc;

        if (isInside()) {
            if (proj > 0.0)
                continue;
        } else {
            if (proj <= 0.0)
                continue;
        }

        stdStringList cc(3);
        cc[0] = toString(cC[n].x());
        cc[1] = toString(cC[n].y());
        cc[2] = toString(cC[n].z());

        std::string::size_type posx = 0;
        while ((posx = of.find_first_of("X", posx)) != std::string::npos) {
            of = of.replace(posx - 1, 3, cc[0]);
        }
        std::string::size_type posy = 0;
        while ((posx = of.find_first_of("Y", posx)) != std::string::npos) {
            of = of.replace(posx - 1, 3, cc[1]);
        }
        std::string::size_type posz = 0;
        while ((posx = of.find_first_of("Z", posx)) != std::string::npos) {
            of = of.replace(posx - 1, 3, cc[2]);
        }
        std::string::size_type posp = 0;
        while ((posp = of.find_first_of("P", posp)) != std::string::npos) {
            of = of.replace(posp - 1, 3, toString(constant::mathConstants::pi));
        }
        formulaParsing realFormula = formulaParsing(of); // f(cc[n](x,y,z))

        real realValue = realFormula.parsing(); // f

        cellQ[n] = realValue;
    }

    return true;
}

bool OpenHurricane::planeRegion::distributing(vectorGeometryArray<cellMesh> &cellQ,
                                              std::string &value) const {
    const auto &mesh = cellQ.mesh();
    std::string of = value;

    if (!isOptionSet()) {
        PLWarning("Option unset for mark region: %d", id());
        return false;
    }

    vector sc = normal_.normalized();
    vector pp = planePoint_;

    const auto &cC = mesh.cellCentre();

    stdStringList orivf;
    of = of.substr(1, of.size() - 2);
    split(of, orivf, ",");

    for (integer n = 0; n < mesh.nCells(); ++n) {
        stdStringList vf = orivf;
        real proj = (cC[n] - pp) * sc;
        if (isInside()) {
            if (proj > 0.0)
                continue;
        } else {
            if (proj <= 0.0)
                continue;
        }

        stdStringList cc(3);
        cc[0] = toString(cC[n].x());
        cc[1] = toString(cC[n].y());
        cc[2] = toString(cC[n].z());

        List<formulaParsing> vecFormula;
        vecFormula.resize(3);
        for (integer i = 0; i < vf.size(); i++) {
            std::string::size_type posx = 0;
            while ((posx = vf[i].find_first_of("X", posx)) != std::string::npos) {
                vf[i] = vf[i].replace(posx - 1, 3, cc[0]);
            }
            std::string::size_type posy = 0;
            while ((posy = vf[i].find_first_of("Y", posy)) != std::string::npos) {
                vf[i] = vf[i].replace(posy - 1, 3, cc[1]);
            }
            std::string::size_type posz = 0;
            while ((posz = vf[i].find_first_of("Z", posz)) != std::string::npos) {
                vf[i] = vf[i].replace(posz - 1, 3, cc[2]);
            }
            std::string::size_type posp = 0;
            while ((posp = vf[i].find_first_of("P", posp)) != std::string::npos) {
                vf[i] = vf[i].replace(posp - 1, 3, toString(constant::mathConstants::pi));
            }
            vecFormula[i] =
                formulaParsing(vf[i]); // x(cc[n](x,y,z))��y(cc[n](x,y,z))��z(cc[n](x,y,z))
        }

        vector vecValue(Zero);
        vecValue.x() = vecFormula[0].parsing();    // u
        vecValue.y() = vecFormula[1].parsing();    // v
        vecValue.z() = vecFormula[2].parsing();    // w

        cellQ[n] = vecValue;
    }

    return true;
}
