/*!
 * \file monitor.hpp
 * \brief Headers of class of monitor.
 *        The subroutines and functions are in the <i>monitor.cpp</i> file.
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
#include "objectFactory.hpp"
#include "runtimeMesh.hpp"

namespace OpenHurricane {
    /*!\brief The base class of monitor.*/
    class monitor {
    private:
        const iteration &iter_;

        const runtimeMesh &mesh_;

    protected:
        /** \brief The step to update monitor. */
        integer updateStep_;

        /** \brief Should print the info. Default is true */
        bool printToScreen_;

        /** \brief Should write the information to the file. Default is off. */
        bool writeToFile_;

        string monitorName_;

        /** \brief Defalut is case name + "_" + monitor name + ".dat" */
        fileName outFile_;

        void removeDuplicate(stringList &residualsNameList) const;

    public:
        declareClassName(monitor);

        declareObjFty(monitor, controller,
                      (const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                       const string &name),
                      (iter, mesh, cont, name));

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        monitor(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                const string &name);

        hur_nodiscard static uniquePtr<monitor> creator(const iteration &iter,
                                                        const runtimeMesh &mesh,
                                                        const controller &cont, const string &name);

        /**
         * \brief Destructor.
         */
        virtual ~monitor() noexcept {}

        hur_nodiscard inline const iteration &iter() const noexcept { return iter_; }

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept { return mesh_; }

        hur_nodiscard inline const string &monitorName() const noexcept { return monitorName_; }

        /** \brief Defalut is case name + "_" + monitor name + ".dat" */
        hur_nodiscard inline const fileName &outFile() const noexcept { return outFile_; }

        /**
         * \brief Monitoring.
         */
        virtual void monitoring() const = 0;
        virtual void subMonitoring() const = 0;
    };
} // namespace OpenHurricane
