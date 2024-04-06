/*!
 * \file relayIN.hpp
 * \brief Header of reading relay files
 *       The subroutines and functions are in the <i>relayIN.cpp</i> file.
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
#pragma once
#include "ArrayInclude.hpp"
#include "dataStructure.hpp"
#include "hdf5I.hpp"
#include "hdf5O.hpp"
#include "relay.hpp"

namespace OpenHurricane {
    class iteration;

    /*!\brief The class of reading relay files.*/
    class relayIN : public relay {
    private:
        iteration &iter_;

        bool isInterpolation_;

        /**
         * \brief If relay file is steady or unsteady.
         * false for steady, true for unsteady
         */
        bool steadyOrUnsteady_;

        using attN = relay::attriNameFile;

    public:
        relayIN(iteration &iter);

        relayIN(const relayIN &) = delete;
        relayIN &operator=(const relayIN &) = delete;

        virtual ~relayIN() noexcept {}

        void readRelay(const fileName &restartFrom);
        void readRelay(const hdf5I &fos, const fileName &restartFrom);

        void createIndexMap(const hdf5I &fos);

        hur_nodiscard inline bool isSteadyRelay() const noexcept {
            return steadyOrUnsteady_ == false;
        }

        hur_nodiscard inline bool isUnsteadyRelay() const noexcept {
            return steadyOrUnsteady_ == true;
        }

    private:
        hur_nodiscard integer getIntegerAttributeFromFile(const hdf5I &fos,
                                                          const string &attrName) const;
        hur_nodiscard std::string getStringAttributeFromFile(const hdf5I &fos,
                                                             const string &attrName) const;
        hur_nodiscard real getRealAttributeFromFile(const hdf5I &fos, const string &attrName) const;

        hur_nodiscard bool fileExist(const hdf5I &fos, const string &attrName) const;

        hur_nodiscard realArray getRealArrayFromFile(const hdf5I &fos,
                                                     const string &attrName) const;
        hur_nodiscard std::string getStringFromFile(const hdf5I &fos, const string &attrName) const;

        void setPhysicalTimeInfo(const real totalPhyTimes, const realArray &lastTimeStep);
        void getPhysicalTimeInfoFromFile(const hdf5I &fos);
        hur_nodiscard stringList getVarName(const hdf5I &fos) const;
    };
} // namespace OpenHurricane
