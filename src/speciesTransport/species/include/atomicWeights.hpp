﻿/*!
 * \file atomicWeights.hpp
 * \brief All the information about the definition of the atomic weights.
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

#include "real.hpp"

namespace OpenHurricane {

    /**
     * \brief The class of the atomic weights.
     */
    class atomicWeight {
    private:
        /**
         * \brief Check tha given atomic is in the table.
         * \return The number of atomic name[3]
         * \retval -1 if not match
         */
        hur_nodiscard static integer checkNum(const char *name);

        hur_nodiscard static bool mapName(const char *tableName, const char *name);

    public:
        /**\brief Structure to hold the element name and atomic weight pair.*/
        struct nameAndWeight {
            char name_[3];
            real weight_;
        };

        static constexpr int nElements_ = 110;

        /**\brief Static table of the weights of all known elements [g/mol or kg/kmol].*/
        static constexpr nameAndWeight atomicWeights_[nElements_] = {
            {"E ", 0},         {"e ", 5.45e-4},   {"H ", 1.00797},   {"D ", 2.01410},
            {"T ", 3.01604},   {"He", 4.00260},   {"Li", 6.93900},   {"Be", 9.01220},
            {"B ", 10.81100},  {"C ", 12.01115},  {"N ", 14.00670},  {"O ", 15.99940},
            {"F ", 18.99840},  {"Ne", 20.18300},  {"Na", 22.98980},  {"Mg", 24.31200},
            {"Al", 26.98150},  {"Si", 28.08600},  {"P ", 30.97380},  {"S ", 32.06400},
            {"Cl", 35.45300},  {"Ar", 39.94800},  {"K ", 39.10200},  {"Ca", 40.08000},
            {"Sc", 44.95600},  {"Ti", 47.90000},  {"V ", 50.94200},  {"Cr", 51.99600},
            {"Mn", 54.93800},  {"Fe", 55.84700},  {"Co", 58.93320},  {"Ni", 58.71000},
            {"Cu", 63.54000},  {"Zn", 65.37000},  {"Ga", 69.72000},  {"Ge", 72.59000},
            {"As", 74.92160},  {"Se", 78.96000},  {"Br", 79.90090},  {"Kr", 83.80000},
            {"Rb", 85.47000},  {"Sr", 87.62000},  {"Y ", 88.90500},  {"Zr", 91.22000},
            {"Nb", 92.90600},  {"Mo", 95.94000},  {"Tc", 99.00000},  {"Ru", 101.07000},
            {"Rh", 102.90500}, {"Pd", 106.40000}, {"Ag", 107.87000}, {"Cd", 112.40000},
            {"In", 114.82000}, {"Sn", 118.69000}, {"Sb", 121.75000}, {"Te", 127.60000},
            {"I ", 126.90440}, {"Xe", 131.30000}, {"Cs", 132.90500}, {"Ba", 137.34000},
            {"La", 138.91000}, {"Ce", 140.12000}, {"Pr", 140.90700}, {"Nd", 144.24000},
            {"Pm", 145.00000}, {"Sm", 150.35000}, {"Eu", 151.96000}, {"Gd", 157.25000},
            {"Tb", 158.92400}, {"Dy", 162.50000}, {"Ho", 164.93000}, {"Er", 167.26000},
            {"Tm", 168.93400}, {"Yb", 173.04000}, {"Lu", 174.99700}, {"Hf", 178.49000},
            {"Ta", 180.94800}, {"W ", 183.85000}, {"Re", 186.20000}, {"Os", 190.20000},
            {"Ir", 192.20000}, {"Pt", 195.09000}, {"Au", 196.96700}, {"Hg", 200.59000},
            {"Tl", 204.37000}, {"Pb", 207.19000}, {"Bi", 208.98000}, {"Po", 210.00000},
            {"At", 210.00000}, {"Rn", 222.00000}, {"Fr", 223.00000}, {"Ra", 226.00000},
            {"Ac", 227.00000}, {"Th", 232.03800}, {"Pa", 231.00000}, {"U ", 238.03000},
            {"Np", 237.00000}, {"Pu", 242.00000}, {"Am", 243.00000}, {"Cm", 247.00000},
            {"Bk", 249.00000}, {"Cf", 251.00000}, {"Es", 254.00000}, {"Fm", 253.00000},
            {"U6", 244.00000}, {"U5", 12.01100},  {"U1", 9.01218},   {"U2", 10.81100},
            {"U3", 24.30500},  {"U4", 26.98154}};

        atomicWeight() {}

        /**
         *\brief Search the atomic weight of atomic name[3]
         * \return The flag of the search operation.
         * \retval - 1 if success,
         * \retval - 0 if false.
         */
        static integer search(const char *name, real &weight);
    };
} // namespace OpenHurricane