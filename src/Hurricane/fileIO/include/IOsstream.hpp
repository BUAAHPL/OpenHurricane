/*!
 * \file IOsstream.hpp
 * \brief Headers of input and output string stream.
 *        The subroutines and functions are in the <i>IOsstream.inl</i> file.
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
#include "errorAbort.hpp"
#include "fileName.hpp"
#include "integer.hpp"
#include <fstream>
#include <iostream>

namespace OpenHurricane {

    /*!
     * \brief The base class of input and output string stream.
     */
    class IOsstream {
    public:
        enum streamFormat { ASCII_FORMAT, BINARY_FORMAT };

        enum accessState { OPENED, CLOSED };

        enum openOptions : short { ONLY_MASTER, CREATE_MASTER, CREATE_ALL };

    private:

        /*! \brief The name of the input file.*/
        fileName name_;

        /*! \brief The stream format.*/
        streamFormat format_;

        /*! \brief The open mode.*/
        std::ios_base::openmode mode_;

        /*! \brief The access state of the file.*/
        accessState openClosed_;

    protected:

        /*! \brief Set file opened.*/
        inline void setOpened();

        /*! \brief Set file opened.*/
        inline void setClosed();

    public:

        /*! \brief Construct null*/
        inline IOsstream();

        /*! \brief Construct setting format and open mode.*/
        inline IOsstream(const streamFormat format, std::ios_base::openmode mode);

        /*! \brief Construct from file name and setting format and open mode.*/
        inline IOsstream(const fileName &name, const streamFormat format = ASCII_FORMAT,
                         std::ios_base::openmode mode = std::ios_base::in);


        /*!\brief Return the const access to the name of the input file.*/
        hur_nodiscard inline const fileName &name() const noexcept;

        /*!\brief Return the non-const access to the name of the input file.*/
       hur_nodiscard inline fileName &name() noexcept;

        /*!\brief Return access to the format.*/
       hur_nodiscard inline streamFormat format() const noexcept;

        /*!\brief Return access to the open mode.*/
       hur_nodiscard inline std::ios_base::openmode mode() const noexcept;

        /*!\brief Return true if the file opened.*/
       hur_nodiscard inline bool opened() const noexcept;

        /*!\brief Return true if the file closed.*/
       hur_nodiscard inline bool closed() const noexcept;

        /*!\brief Change the open mode.*/
       inline void changeMode(const std::ios_base::openmode _mode) noexcept;

        /*!\brief Change the open file name.*/
        inline void changeFileName(const fileName &_name) noexcept;

        /*!\brief Change the format.*/
        inline void changeFormat(const streamFormat _format) noexcept;

        /*! \brief Return true if the file is already existed.*/
        hur_nodiscard inline bool existed() const;
    };

} // namespace OpenHurricane

#include "IOsstream.inl"