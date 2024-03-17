/*!
 * \file dataSections.hpp
 * \brief Header of data sections
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
#pragma once

namespace OpenHurricane {
namespace dataSections {
extern const char magicNumber[];

extern const int unitIntFlag;

extern const int integerFlag;

extern const int ndFlag;

extern const float gridSection;

extern const float nodeTotal;
extern const float nodeZone;

extern const float cellTotal;
extern const float cellZone;

extern const float faceTotal;
extern const float faceZone;

extern const float dataHeader;
extern const float dataSection;
} // namespace dataSections

namespace Options {
namespace writeRelayOption {
enum writeRelayOptions : short {
    ONLY_DATA = 0,   /*!< \brief Only write data to relay file. */
    DATA_GRID = 1,   /*!< \brief Write data and grid to relay file. */
    ONLY_GRID = 2,   /*!< \brief Only write grid to relay file. */
    CASE_CONFIG = 3, /*!< \brief Only write grid and controller to case file. */
    NO_WRITE = 4     /*!< \brief No write. */
};
}
} // namespace Options
} // namespace OpenHurricane
