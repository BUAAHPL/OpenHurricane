/*!
 * \file processTopology.hpp
 * \brief Headers of the process topology.
 *        The subroutines and functions are in the <i>processTopology.cpp</i> file.
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
#include "zoneMesh.hpp"

namespace OpenHurricane {
    template <class zoneList> class processTopology {
    private:
        const zoneList &zones_;

        /**
         * \brief The index of zone that send data from this process.
         */
        integerListList zoneIdSend_;

        /**
         * \brief The index of zone that receive data on this process from the remote process.
         */
        integerListList zoneIdRecv_;

        /**
         * \brief Get zone id.
         * \param[in] nLayer - The number of ghost cell layers
         */
        void getZoneId(const integer nLayer);

    public:

        /**
         * \brief Disallow null constructor.
         */
        processTopology() = delete;

        /**
         * \brief Construct from zone list and the number of ghost cell layers.
         */
        processTopology(const zoneList &zones, const integer nLayer);

        /**
         * \brief Disallow copy constructor.
         */
        processTopology(const processTopology &) = delete;

        /**
         * \brief Destructor.
         */
        inline ~processTopology() noexcept;

        /**
         * \brief Zone list.
         */
        hur_nodiscard inline const zoneList &zones() const noexcept;

        /**
         * \brief The index of zone that send data from this process.
         */
        hur_nodiscard inline const integerListList &zoneIdSend() const noexcept;

        /**
         * \brief The index of zone that receive data on this process from the remote process.
         */
        hur_nodiscard inline const integerListList &zoneIdRecv() const noexcept;

        hur_nodiscard inline integer zoneIdSend(const integer nLayer, const integer n) const;
        hur_nodiscard inline integer zoneIdRecv(const integer nLayer, const integer n) const;

        /**
         * \brief Disallow bitwise assignment.
         */
        void operator=(const processTopology &) = delete;
    };

    
    template <class zoneType, class dataType> class processSend {
    private:
        const zoneType &zone_;

        Array<typename feature<dataType>::elementType> buf_;

    public:
        processSend() = delete;

        inline processSend(const zoneType &zone);

        processSend(const processSend &) = delete;

        inline ~processSend() noexcept;

        inline void clear() noexcept;

        void writeBUf(const Array<dataType> &f);

        void send() const;

        void sendNonBlock(HurMPI::Request *request);

        void operator=(const processSend &) = delete;
    };
        
    template <class zoneType, class dataType> class processReceiv {
    private:
        const zoneType &zone_;

        Array<typename feature<dataType>::elementType> buf_;

    public:
        processReceiv() = delete;

        inline processReceiv(const zoneType &zone);

        processReceiv(const processReceiv &) = delete;

        inline ~processReceiv() noexcept;

        inline void clear() noexcept;

        void readBUf(Array<dataType> &f);

        void recv();

        hur_nodiscard inline const zoneType &zone() const noexcept { return zone_; }

        void recvNonBlock(HurMPI::Request *request);
        void recvNonBlock(HurMPI::Request *request, const integer sizeBuf);

        hur_nodiscard inline Array<typename feature<dataType>::elementType> &buf() noexcept;

        void operator=(const processReceiv &) = delete;
    };
} // namespace OpenHurricane

#include "processTopology.inl"