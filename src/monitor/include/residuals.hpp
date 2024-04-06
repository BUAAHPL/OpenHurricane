/*!
 * \file residuals.hpp
 * \brief Headers of class of monitoring residuals.
 *        The subroutines and functions are in the <i>residuals.cpp</i> file.
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
#include "monitor.hpp"

namespace OpenHurricane {
    /*!\brief The class of monitoring residuals.*/
    class residuals : public monitor {
    private:
        /** \brief The name list of variables monitored. */
        stringList residualsNameList_;

        stringList residualsCmptNameList_;

        /** \brief Check convergence by monitoring residuals? */
        bool checkResidualConvergence_;

        /** \brief The convergence threshold for each variables monitored. */
        realArray convergenceThreshold_;

        integerList resCmptNum_;

        /**
         * \brief Read from controller.
         */
        void readFromCont(const controller &cont);

        void setDefaultMonitorVars();

        mutable fileOsstream fosResidual_;

        void parsingList(const controller &resCont, const integerList &cmptNum,
                         stringList &cmptName);

    public:
        declareClassName(residuals);

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        residuals(const iteration &iter, const runtimeMesh &mesh, const controller &cont,
                  const string &name);

        /**
         * \brief Destructor.
         */
        virtual ~residuals() noexcept { fosResidual_.close(); }

        bool resConvergence() const;

        void writeResiduals(fileOsstream &fos) const;

        void writeSubIterResiduals(fileOsstream &fos) const;

        /** \brief The name list of variables monitored. */
        hur_nodiscard inline const stringList &residualsNameList() const noexcept {
            return residualsNameList_;
        }

        /** \brief Check convergence by monitoring residuals? */
        hur_nodiscard inline bool checkResidualConvergence() const noexcept {
            return checkResidualConvergence_;
        }

        /** \brief The convergence threshold for each variables monitored. */
        hur_nodiscard inline const realArray &convergenceThreshold() const noexcept {
            return convergenceThreshold_;
        }

        /**
         * \brief Monitoring.
         */
        inline virtual void monitoring() const {
            if (iter().cStep() % updateStep_ == 0) {
                writeResiduals(fosResidual_);
            }
        }
        inline virtual void subMonitoring() const {
            if (iter().subIter().cSubStep() % updateStep_ == 0) {
                writeSubIterResiduals(fosResidual_);
            }
        }
    };
} // namespace OpenHurricane
