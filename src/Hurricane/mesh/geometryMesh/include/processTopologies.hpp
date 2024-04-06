/*!
 * \file processTopologies.hpp
 * \brief Headers of the process topology.
 *        The subroutines and functions are in the <i>processTopologies.cpp</i> file.
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

#include "processTopology.hpp"

namespace OpenHurricane {
    using processCutTopology = processTopology<cutZoneList>;
    using processPerTopology = processTopology<perZoneList>;

    template <class dataType> using cutProcessSend = processSend<cutZone, dataType>;
    template <class dataType> using perProcessSend = processSend<perZone, dataType>;

    template <> void cutProcessSend<real>::writeBUf(const Array<real> &f);
    template <> void perProcessSend<real>::writeBUf(const Array<real> &f);

    template <> void cutProcessSend<integer>::writeBUf(const Array<integer> &f);
    template <> void perProcessSend<integer>::writeBUf(const Array<integer> &f);

    template <> void cutProcessSend<complex>::writeBUf(const Array<complex> &f);
    template <> void perProcessSend<complex>::writeBUf(const Array<complex> &f);

    template <> void cutProcessSend<realArray>::writeBUf(const Array<realArray> &f);
    template <> void perProcessSend<realArray>::writeBUf(const Array<realArray> &f);

    template <class dataType> class periodicTransferSameProc {
    private:
        const perZone &zone_;

    public:
        periodicTransferSameProc() = delete;

        inline periodicTransferSameProc(const perZone &zone) : zone_(zone) {
            if (!zone_.isSameProc()) {
                LFatal("Periodic zone must be in the same process");
            }
        }

        periodicTransferSameProc(const periodicTransferSameProc &other) = delete;
        periodicTransferSameProc &operator=(const periodicTransferSameProc &other) = delete;

        inline ~periodicTransferSameProc() noexcept {}

        inline void transferData(const Array<dataType> &var, Array<dataType> &recvVar) {
            const auto &sorId = zone_.sor();
            const auto &desId = zone_.des();

            for (integer i = 0; i < sorId.size(); ++i) {
                recvVar[desId[i]] = var[sorId[i]];
            }
        }
    };

    template <>
    inline void periodicTransferSameProc<vector>::transferData(const Array<vector> &var,
                                                               Array<vector> &recvVar) {
        const auto &sorId = zone_.sor();
        const auto &desId = zone_.des();

        for (integer i = 0; i < sorId.size(); ++i) {
            if (zone_.isRotational()) {
                recvVar[desId[i]] = zone_.RotionalMatrix() * var[sorId[i]];
            } else {
                recvVar[desId[i]] = var[sorId[i]];
            }
        }
    }

    template <class dataType> using cutProcessReceiv = processReceiv<cutZone, dataType>;

    template <class dataType> using perProcessReceiv = processReceiv<perZone, dataType>;

    template <> void cutProcessReceiv<real>::readBUf(Array<real> &f);
    template <> void perProcessReceiv<real>::readBUf(Array<real> &f);
    template <> void cutProcessReceiv<real>::recvNonBlock(HurMPI::Request *request);
    template <> void perProcessReceiv<real>::recvNonBlock(HurMPI::Request *request);

    template <> void cutProcessReceiv<integer>::readBUf(Array<integer> &f);
    template <> void perProcessReceiv<integer>::readBUf(Array<integer> &f);
    template <> void cutProcessReceiv<integer>::recvNonBlock(HurMPI::Request *request);
    template <> void perProcessReceiv<integer>::recvNonBlock(HurMPI::Request *request);

    template <> void cutProcessReceiv<complex>::readBUf(Array<complex> &f);
    template <> void perProcessReceiv<complex>::readBUf(Array<complex> &f);
    template <> void cutProcessReceiv<complex>::recvNonBlock(HurMPI::Request *request);
    template <> void perProcessReceiv<complex>::recvNonBlock(HurMPI::Request *request);

    template <> void cutProcessReceiv<realArray>::readBUf(Array<realArray> &f);
    template <> void perProcessReceiv<realArray>::readBUf(Array<realArray> &f);
    template <> void cutProcessReceiv<realArray>::recvNonBlock(HurMPI::Request *request);
    template <> void perProcessReceiv<realArray>::recvNonBlock(HurMPI::Request *request);
    template <>
    void cutProcessReceiv<realArray>::recvNonBlock(HurMPI::Request *request, const integer sizeBuf);
    template <>
    void perProcessReceiv<realArray>::recvNonBlock(HurMPI::Request *request, const integer sizeBuf);

    template <> void perProcessReceiv<vector>::readBUf(Array<vector> &f);

    template <> void perProcessReceiv<tensor>::readBUf(Array<tensor> &f);
} // namespace OpenHurricane
