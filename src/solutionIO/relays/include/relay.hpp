/*!
 * \file relay.hpp
 * \brief Header of relay files
 *       The subroutines and functions are in the <i>relay.cpp</i> file.
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

namespace OpenHurricane {
    namespace relayDefine {
        enum writeRelayOptions {
            ONLY_DATA = 0,   /*!< \brief Only write data to relay file. */
            DATA_GRID = 1,   /*!< \brief Write data and grid to relay file. */
            ONLY_GRID = 2,   /*!< \brief Only write grid to relay file. */
            CASE_CONFIG = 3, /*!< \brief Only write grid and controller to case file. */
            NO_WRITE = 4     /*!< \brief No write. */
        };

        enum state { steady = 0, unsteady = 1 };
    } // namespace relayDefine

    class iteration;

    /*!\brief The base class of relay files.*/
    class relay {
    public:
        /**
         * \brief Attribute name for or from file.
         */
        class attriNameFile {
        public:
            static const string program;
            static const string version;
            static const string dateTime;
            static const string meshFile;
            static const string nProcessor;
            static const string dataType;
            static const string totalStep;
            static const string state;
            static const string time;
            static const string lastTimeStep;
            static const string nVariables;
            static const string variableName;
            static const string cellOriginIndex;
            static const string cellCentre;
            static const string nTimeGroups;
            static const string timeGroup0;
            hur_nodiscard inline static string timeGroupm(const integer i) {
                return string("timeGroupm") + toString(i);
            }
        };

    public:
        inline relay(){};

        relay(const relay &) = delete;
        relay &operator=(const relay &) = delete;

        virtual ~relay() noexcept {}
    };
} // namespace OpenHurricane
