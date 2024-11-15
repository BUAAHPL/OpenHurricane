﻿/*!
 * \file relayOUT.hpp
 * \brief Header of writting relay files
 *       The subroutines and functions are in the <i>relayOUT.cpp</i> file.
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

    /*!\brief The class of writting relay files.*/
    class relayOUT : public relay {
    private:
        const iteration &iter_;
        mutable hdf5O ofos_;

        using attN = relay::attriNameFile;

    public:
        relayOUT(const iteration &iter);

        relayOUT(const relayOUT &) = delete;
        relayOUT &operator=(const relayOUT &) = delete;

        virtual ~relayOUT() noexcept {}

        void writeRelay(const fileName &relayFileName, const fileName &meshFileName) const;

        /*!\brief Write object's original index to relay file.*/
        void writeOriIndex(hdf5O &fos) const;

        void writeCellCentre(hdf5O &fos) const;

    private:
        inline void writeStringAttriToFile(const std::string &strAtt,
                                           const string &attrName) const {
            if (HurMPIBase::master()) {
                ofos_.writeStringAttributeToFile(strAtt, attrName);
            }
        }
        inline void writeIntAttriToFile(const integer intAtt, const string &attrName) const {
            if (HurMPIBase::master()) {
                ofos_.writeIntegerAttributeToFile(intAtt, attrName);
            }
        }
    };
} // namespace OpenHurricane
