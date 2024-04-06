/*!
 * \file shockSensor.cpp
 * \brief Main subroutines for shock sensor.
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

#include "shockSensor.hpp"

OpenHurricane::shockSensor::shockSensor(const flowModel &flows, const controller &cont)
    : flows_(flows), delta_(flows.mesh().nCells(), Zero),
      sensor_(object("shockSensor", flows.mesh(), object::NOT_WRITE), flows.mesh()),
      thetaIJ_(flows.mesh().nFaces(), Zero),
      differThr_(cont.findOrDefault<real>("perssureSensorK", 1.0 / 128.0)),
      differTThr_(cont.findOrDefault<real>("temperatureSensorK", 0.002)), checkTemperature_(false),
      THighLimit_(cont.findOrDefault<real>("THighLimitShock", 100.0)) {
    const integer nCells = flows_.mesh().nCells();

    for (integer n = 0; n < nCells; ++n) {
        delta_[n] = pow(flows_.mesh().cellVol()[n], real(1.0 / 3.0));
    }
    sensor_ = Zero;
    getGeometricalWeights();

    controllerSwitch myContS(cont);
}

void OpenHurricane::shockSensor::calcShockSensor() {
    const integer nCells = flows_.mesh().nCells();
    const auto &cv = flows_.mesh().cellVolume();
    sensor_ = Zero;
    for (integer n = 0; n < nCells; ++n) {
        const real S = skewMagnitude(v().grad()[n]);

        const real dil = v().grad()[n].xx() + v().grad()[n].yy() + v().grad()[n].zz();

        const real c = sqrt(flows_.gama()[n] * flows_.p()[n] / flows_.rho()[n]);
        const real dil2 = sqr(dil);
        const real CSS =
            0.5 *
            (1.0 - tanh(real(2.5 + 10.0 * pow(cv[n], real(1.0 / 3.0)) / max(c, tiny) * dil))) *
            dil2 / (dil2 + sqr(S) + tiny);

        if (CSS > real(0.02)) {
            sensor_[n] = 1.0;
        }
    }

    if (checkTemperature_) {
        checkTemperatureDiffer();
    }
    const auto &fZ = flows_.mesh().faceZones();
    const auto &fL = flows_.mesh().faces();
    realTransfer myTransfer(flows_.mesh(), sensor_, false, true);
    myTransfer.transferInit();
    for (integer fZI = 0; fZI < fZ.size(); ++fZI) {
        if (fZ[fZI].isInterior()) {
            continue;
        } else if (fZ[fZI].isCutFace()) {
            continue;
        } else if (fZ[fZI].isPeriodic() || fZ[fZI].isPeriodicShadow()) {
            continue;
        } else {
            for (integer fI = fZ[fZI].firstIndex(); fI <= fZ[fZI].lastIndex(); ++fI) {
                const auto &cl = fL[fI].leftCell();
                const auto &cr = fL[fI].rightCell();
                sensor_[cr] = sensor_[cl];
            }
        }
    }
    myTransfer.transferring();
}

void OpenHurricane::shockSensor::getGeometricalWeights() {
    const auto &mesh = flows_.mesh();
    const auto &fl = mesh.faces();
    const auto &cells = mesh.cells();
    const auto &cC = mesh.cellCentre();

    for (integer i = 0; i < mesh.nCells(); ++i) {
        const auto &ci = cC[i];
        tensor I = Zero;
        vector Rs = Zero;

        for (integer j = 0; j < cells[i].faceSize(); ++j) {
            const auto fi = cells[i].facei(j);
            const auto cl = fl[fi].leftCell();
            const auto cr = fl[fi].rightCell();
            const integer jj = cl + cr - i;
            vector tempV = (cC[jj] - ci);
            I += inv(tempV.magSqr()) * tempV & tempV;
            Rs += tempV;
        }
        vector lamb = Rs / I;

        for (integer j = 0; j < cells[i].faceSize(); ++j) {
            const auto fi = cells[i].facei(j);
            const auto cl = fl[fi].leftCell();
            const auto cr = fl[fi].rightCell();
            const integer jj = cl + cr - i;
            vector tempV = (cC[jj] - ci);
            real thet = real(1.0) + inv(tempV.magSqr()) * lamb * tempV;
            thet = min(real(2), max(thet, real(0)));

            if (i == cl) {
                thetaIJ_[fi][0] = thet;
            } else {
                thetaIJ_[fi][1] = thet;
            }
        }
    }
}

void OpenHurricane::shockSensor::checkDiffer(const realArray &p, const real differThr) {
    const auto &mesh = flows_.mesh();

    const integer nCells = mesh.nCells();
    const auto &faces = mesh.faces();
    const auto &cells = mesh.cells();

    for (integer celli = 0; celli < nCells; ++celli) {
        if (sensor_[celli] == 1.0) {
            continue;
        }
        // numerator
        real num = Zero;

        // denominator
        real den = Zero;

        for (integer i = 0; i < cells[celli].faceSize(); ++i) {
            const integer fi = cells[celli].facei(i);
            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();

            const integer I = celli;
            const integer J = (cl == celli) ? cr : cl;
            if (cl == celli) {
                num += (thetaIJ_[fi][0] * (p[J] - p[I]));
            } else {
                num += (thetaIJ_[fi][1] * (p[J] - p[I]));
            }
            den += (p[J] + p[I]);
        }
        const real pgI = 0.5 * mag(num) / max(den, tiny);

        if (pgI >= differThr) {
            sensor_[celli] = 1.0;
        } else {
            sensor_[celli] = pgI / max(differThr, tiny);
            if (sensor_[celli] > 1) {
                sensor_[celli] = 1;
            }
        }
    }
}

void OpenHurricane::shockSensor::checkTemperatureDiffer() {
    const auto &mesh = flows_.mesh();

    const integer nCells = mesh.nCells();
    const auto &faces = mesh.faces();
    const auto &cells = mesh.cells();

    for (integer celli = 0; celli < nCells; ++celli) {
        if (sensor_[celli] == 1.0 || flows_.T()[celli] < THighLimit_) {
            continue;
        }
        // numerator
        real num = Zero;

        // denominator
        real den = Zero;

        for (integer i = 0; i < cells[celli].faceSize(); ++i) {
            const integer fi = cells[celli].facei(i);
            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();

            const integer I = celli;
            const integer J = (cl == celli) ? cr : cl;
            if (cl == celli) {
                num += (thetaIJ_[fi][0] * (flows_.T()[J] - flows_.T()[I]));
            } else {
                num += (thetaIJ_[fi][1] * (flows_.T()[J] - flows_.T()[I]));
            }
            den += (flows_.T()[J] + flows_.T()[I]);
        }
        const real pgI = 0.5 * mag(num) / max(den, tiny);

        if (pgI >= differTThr_) {
            sensor_[celli] = 1.0;
        } else {
            /*sensor_[celli] = pgI / max(differTThr_, tiny);
            if (sensor_[celli] > 1)
            {
                sensor_[celli] = 1;
            }*/
        }
    }
}

void OpenHurricane::shockSensor::checkTemperatureDiffer(realArray &factor) const {
    const auto &mesh = flows_.mesh();

    const integer nCells = mesh.nCells();
    const auto &faces = mesh.faces();
    const auto &cells = mesh.cells();

    for (integer celli = 0; celli < nCells; ++celli) {
        if (sensor_[celli] == 1.0 || flows_.T()[celli] < THighLimit_) {
            continue;
        }
        // numerator
        real num = Zero;

        // denominator
        real den = Zero;

        for (integer i = 0; i < cells[celli].faceSize(); ++i) {
            const integer fi = cells[celli].facei(i);
            const auto cl = faces[fi].leftCell();
            const auto cr = faces[fi].rightCell();

            const integer I = celli;
            const integer J = (cl == celli) ? cr : cl;
            if (cl == celli) {
                num += (thetaIJ_[fi][0] * (flows_.T()[J] - flows_.T()[I]));
            } else {
                num += (thetaIJ_[fi][1] * (flows_.T()[J] - flows_.T()[I]));
            }
            den += (flows_.T()[J] + flows_.T()[I]);
        }
        const real pgI = 0.5 * mag(num) / max(den, tiny);

        if (pgI >= differTThr_) {
            factor[celli] = 1.0;
        } else {
            factor[celli] = max(factor[celli], pgI / max(differTThr_, tiny));
            if (factor[celli] > 1) {
                factor[celli] = 1;
            }
        }
    }
}