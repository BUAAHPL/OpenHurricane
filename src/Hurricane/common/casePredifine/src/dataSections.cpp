/*!
 * \file dataSections.cpp
 * \brief The subroutines and functions of data sections
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

#include "dataSections.hpp"

const char OpenHurricane::dataSections::magicNumber[] = "#!UNS100";

const int OpenHurricane::dataSections::unitIntFlag = 1;
const int OpenHurricane::dataSections::integerFlag = 101;
const int OpenHurricane::dataSections::ndFlag = 3;

const float OpenHurricane::dataSections::gridSection = 187.0f;

const float OpenHurricane::dataSections::nodeTotal = 100.0f;
const float OpenHurricane::dataSections::nodeZone = 101.0f;

const float OpenHurricane::dataSections::cellTotal = 120.0f;
const float OpenHurricane::dataSections::cellZone = 121.0f;

const float OpenHurricane::dataSections::faceTotal = 130.0f;
const float OpenHurricane::dataSections::faceZone = 131.0f;

const float OpenHurricane::dataSections::dataHeader = 299.0f;
const float OpenHurricane::dataSections::dataSection = 357.0f;