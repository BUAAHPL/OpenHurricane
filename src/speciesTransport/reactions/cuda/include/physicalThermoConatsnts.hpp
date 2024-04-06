/*!
 * \file physicalThermoConatsnts.hpp
 * \brief Header of physical thermo conatsnts in CUDA platform.
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

#ifdef CUDA_PARALLEL
#include "CUDAFunctions.hpp"

#ifndef CUDA_Patm
/**\brief Standard atmosphere pressure (Unit: [Pa]).*/
#define CUDA_Patm cu_real(1.01325e5)
#endif // !CUDA_Patm

#ifndef CUDA_Ru
/**\brief Universal gas constant (Unit: [J/(kmol K)]).*/
#define CUDA_Ru cu_real(8.3144626181532403e3)
#endif // !CUDA_Ru

#endif // CUDA_PARALLEL
