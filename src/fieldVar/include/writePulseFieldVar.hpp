/*!
 * \file writePulseFieldVar.hpp
 * \brief Headers of class of writing pulse field variables.
 *        The subroutines and functions are in the <i>writePulseFieldVar.cpp</i> file.
 * \author Chen Zhenyi
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

#pragma once
#include "timeSumField.hpp"
#include "writeFieldVar.hpp"

namespace OpenHurricane {
    class calculateFieldVar;
    /*!\brief The class of write pulse field variables.*/
    class writePulseFieldVar : public writeFieldVar, public timeSumField {
    private:
        real startTime_;
        real curretTimeGap_;
        real endTimeAveCal_;

        void changeOutVarNameToPulse();

    public:
        declareClassNames;

        writePulseFieldVar(const flowModel &flows, const iteration &iter, const controller &cont,
                           const string &writeId, std::map<std::string, object *> &outFieldVarMap);

        /*!\brief Destructor.*/
        virtual ~writePulseFieldVar() noexcept {}

        virtual void updating();

        virtual void correct() const;
        virtual void reseting() const;

    protected:
        const realArray calcTotalPressure(const string &Ma, const realArray &P,
                                          const realArray &gama, const vectorArray &v,
                                          const realArray &rho);

        const realArray calcTotalTemperature(const string &Ma, const realArray &T,
                                             const realArray &gama, const vectorArray &v,
                                             const realArray &P, const realArray &rho);

    protected:
        template <class Type, class GeometryMesh>
        Array<Type> calcTimeAveVar(const geometryArray<Type, GeometryMesh> &var) const {
            return writeFieldVar::calcTimeAveVar(var, curretTimeGap_);
        }

        void calcTimeAveSpecies();

        void setCurretTimeGap();

        virtual fileName getFileName() const;
        virtual bool writeNow() const;
    };
} // namespace OpenHurricane
