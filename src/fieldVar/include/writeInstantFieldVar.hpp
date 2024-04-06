/*!
 * \file writeInstantFieldVar.hpp
 * \brief Headers of class of writing instant field variables.
 *        The subroutines and functions are in the <i>writeInstantFieldVar.cpp</i> file.
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
#include "writeFieldVar.hpp"

namespace OpenHurricane {
    class calculateFieldVar;
    /*!\brief The class of write instant field variables.*/
    class writeInstantFieldVar : public writeFieldVar {
    public:
        declareClassNames;

        writeInstantFieldVar(const flowModel &flows, const iteration &iter, const controller &cont,
                             const string &writeId,
                             std::map<std::string, object *> &outFieldVarMap);

        /*!\brief Destructor.*/
        virtual ~writeInstantFieldVar() noexcept {}

        virtual void updating();

    protected:
        virtual fileName getFileName() const;
        virtual bool writeNow() const;

        void calcStagnationParameters() const;
        void calcViscousRatio() const;
        void calcVorticity() const;
        void calcDeltaVorticity() const;
        void calcOtherVorticity() const;
        void calcMoleSpecies() const;
    };
} // namespace OpenHurricane
