/*!
 * \file fileIsstream.hpp
 * \brief Headers of file input string stream.
 *        The subroutines and functions are in the <i>fileIsstream.cpp</i> file.
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
#include <fstream>

namespace OpenHurricane {

    /*!
     * \brief The class of file input string stream.
     */
    class fileIsstream : public IOsstream {
    private:
        /*! \brief The buffer string of the file.*/
        mutable uniquePtr<std::string> bufferPtr_;

        /*! \brief The file input stream.*/
        std::ifstream if_;

        /*! \brief Is the file buffered?*/
        bool isBuffered_;

    public:
        /*! \brief Construct null.*/
        fileIsstream();

        /*! \brief Construct setting format and open mode.*/
        fileIsstream(const streamFormat format, std::ios_base::openmode mode);

        /*!
         * \brief Construct from file name and setting format and open mode.
         *        And open the file.
         */
        fileIsstream(const fileName &name, const streamFormat format = ASCII_FORMAT,
                     std::ios_base::openmode mode = std::ios_base::in);

        inline ~fileIsstream() noexcept {
            if (if_.is_open()) {
                close();
            }
            clear();
        }

        /*!\brief Return the const access to the buffer string.*/
        hur_nodiscard inline const std::string &bufferString() const noexcept;

        /*!
         * \brief Return the non-const access to the buffer string.
         *        Use with care.
         */
        hur_nodiscard inline std::string &bufferString();

        /*!\brief Is the file opened.*/
        hur_nodiscard inline bool isOpen() const noexcept;

        /*!\brief Is the file buffered.*/
        hur_nodiscard inline bool isBuffered() const noexcept;

        /*!\brief Return the flag of the file input stream.*/
        hur_nodiscard inline int flags() const noexcept;

        // Open and close file

        void open();

        void close();

        inline void clear() noexcept {
            bufferPtr_.clear();
        }

        /*!\brief Read the whole file to the buffer string*/
        void readBuffer();
    };

} // namespace OpenHurricane

#include "fileIsstream.inl"