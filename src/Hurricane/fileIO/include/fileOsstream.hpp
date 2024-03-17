/*!
 * \file fileOsstream.hpp
 * \brief Headers of file output string stream.
 *        The subroutines and functions are in the <i>fileOsstream.cpp</i> file.
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
#include "IOsstream.hpp"
#include "fileName.hpp"
#include "real.hpp"
#include <fstream>

namespace OpenHurricane {

    /*!
     * \brief The class of file input string stream.
     */
    class fileOsstream : public IOsstream {
    private:
        /*! \brief The file input stream.*/
        std::ofstream of_;

        std::streamsize defaultCoutPrecision_;

        openOptions openOption_;

    public:
        /*! \brief Construct null.*/
        fileOsstream();

        /*! \brief Construct setting format and open mode.*/
        fileOsstream(const streamFormat format, std::ios_base::openmode mode);

        /*!
         * \brief Construct from file name and setting format and open mode.
         *        And open the file.
         */
        fileOsstream(const fileName &name, const streamFormat format = ASCII_FORMAT,
                     std::ios_base::openmode mode = std::ios_base::out);

        /*! \brief Construct setting format and open mode.*/
        fileOsstream(const streamFormat format, const openOptions _op,
                     std::ios_base::openmode mode);

        /*!
         * \brief Construct from file name and setting format and open mode.
         *        And open the file.
         */
        fileOsstream(const fileName &name, const openOptions _op,
                     const streamFormat format = ASCII_FORMAT,
                     std::ios_base::openmode mode = std::ios_base::out);

        /*!\brief Destructor.*/
        inline ~fileOsstream() noexcept {
            if (of_.is_open()) {
                close();
            }
        }

        /*!\brief Is the file opened.*/
        hur_nodiscard inline bool isOpen() const noexcept;

        /*!\brief Return the flag of the file output stream.*/
        hur_nodiscard inline int flags() const noexcept;

        hur_nodiscard inline openOptions openOption() const noexcept;

        // Open and close file

        void open();

        void close();

        std::ofstream &flush();

        hur_nodiscard std::ofstream &os() noexcept;

        std::ofstream &write(const char *_str, std::streamsize _count);

        std::ofstream &write(const std::string &_str, std::streamsize _count);

        std::ofstream &setRealPrecision(const int _precision = feature<real>::precision);
        std::ofstream &unsetRealPrecision();

        std::ofstream &setScientic(const int _precision = feature<real>::precision);
        std::ofstream &unsetScientic();
    };

} // namespace OpenHurricane

#include "fileOsstream.inl"