/*!
 * \file writeFaceZone.hpp
 * \brief Headers of base class of writeFaceZone.
 *        The subroutines and functions are in the <i>writeFaceZone.cpp</i> file.
 * \author Chen Zhenyi
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
#include "calculateFaceZoneVar.hpp"
#include "flowModel.hpp"
#include <set>

namespace OpenHurricane {
    //class calculateFieldVar;
    /*!\brief The class of writeFaceZone.*/
    class writeFaceZone {
    private:
        const flowModel &flows_;
        const controller &cont_;

        stringList writeFaceZoneNameList_;
        List<stringList> writeFaceZoneVarList_;
        integerList writeFaceZoneIdList_;

        void checkDuplicate(const faceZone &fz, const stringList &outV) const;

        void setFaceZoneList();
        void setFaceZoneId();

        void writeToFile() const;

        /**
         * \brief Write to tecplot file.
         */
        void writeToTecplot(fileOsstream &fos, const integer fzid,
                            const stringList &outVarName) const;

        void updating() const;

        fileName getFileName(const string &zone) const;

        /**\brief Map for output variables of flow field.*/
        mutable std::map<std::string, object *> outFieldVarMap_;

        /**\brief Map for output variables of flow field.*/
        mutable std::map<std::string, object *> outSamplingFieldVarMap_;

        /** \brief faceZoneVarMap_[faceZoneName][variableName] */
        std::map<std::string, std::set<std::string>> faceZoneVarMap_;

        void setFaceZoneVarMap();

        void clearOutFieldVarMap() const;

        void resetOutFieldVarMap() const;

        faceZoneFieldParameterFuncMap faceZoneVarFuncMap_;

        bool isSampling_;

        integer samplingStep_;

        /**
         * \brief The time for start output average or pulse solutions. Unit: [s].
         */
        real startTime_;

        mutable real timeElasped_;
        mutable real lastTime_;
        mutable bool firstCalling_;

    public:
        /**
         * \brief Construct from controller and runtimeMesh.
         */
        writeFaceZone(const flowModel &flows, const controller &cont);

        /**
         * \brief Destructor.
         */
        inline ~writeFaceZone() noexcept { clearOutFieldVarMap(); }

        inline const stringList &writeFaceZoneNameList() const noexcept;
        inline const List<stringList> &writeFaceZoneVarList() const noexcept;

        void write() const;

    protected:
        bool samplingNow() const;

        void sampling() const;

        void samplingFromField(const faceZone &fz, const string &varN, const real dt,
                               const integer fzid) const;
        void samplingFromMap(const faceZone &fz, const string &varN, const real dt,
                             const integer fzid) const;

    public:
        bool writeNow() const;
        /**
         * \brief Find the varibale if in the output table.
         * \param[in] fzName - The name of face zone.
         * \param[in] varName - The variable name.
         * \return True if found.
         * \retval A bool value.
         */
        bool findInMap(const string &fzName, const string &varName) const;

        void setOutputField(const faceZone &fz, const string &varName,
                            const realArray &value) const;
        void setOutputField(const faceZone &fz, const string &varName, realArray &&value) const;
        void setOutputField(const faceZone &fz, const string &varName,
                            const vectorArray &value) const;
        void setOutputField(const faceZone &fz, const string &varName, vectorArray &&value) const;

        hur_nodiscard inline faceZoneFieldParameterFuncMap &faceZoneVarFuncMap() noexcept;
        hur_nodiscard inline const faceZoneFieldParameterFuncMap &
        faceZoneVarFuncMap() const noexcept;

    };
} // namespace OpenHurricane
#include "writeFaceZone.inl"