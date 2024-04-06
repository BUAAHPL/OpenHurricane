/*!
 * \file solutionWrite.hpp
 * \brief Headers of class of solution write-out.
 *        The subroutines and functions are in the <i>solutionWrite.cpp</i> file.
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

#include "flowModel.hpp"
#include "runtimeMesh.hpp"
#include "solutionWrites.hpp"
#include <set>

namespace OpenHurricane {
    /*!\brief The base class of solutionWrite.*/
    class solutionWrite {
    public:
        enum class solutionType : short {
            onlyInstantaneous = 0,
            onlyAverage = 1,
            onlyPulse = 2,
            InatsntAndAve = 3,
            InstantAndAveAndPulse = 4
        };

    private:
        const flowModel &flows_;

        const iteration &iter_;

        const runtimeMesh &mesh_;

        std::set<std::string> sets_;

        /**\brief Map for output variables of flow field.*/
        mutable std::map<std::string, object *> outFieldVarMap_;

        mutable std::map<std::string, object *> outAveFieldVarMap_;
        mutable std::map<std::string, object *> outPulseFieldVarMap_;

        void getOutSetMap(const stringList &outFieldVarList);

        solutionType solType_;

        /**
         * \brief The time for start output average or pulse solutions. Unit: [s].
         */
        real startTime_;

        mutable real timeElasped_;
        mutable real lastTime_;
        mutable bool firstCalling_;

        stringList outVarNameList_;

        bool isWriteProfile_;

        bool isSampling_;

        integer samplingStep_;

    protected:
        hur_nodiscard stringList outVarNameList() const;

        void removeDuplicate(stdStringList &writeList) const;

        void getVarList(const controller &cont);

        hur_nodiscard fileName outFile(const char *c = ".dat") const;

        void clearOutVarMap(std::map<std::string, object *> &outFieldVarMap) const;

        void calcAveVar() const;
        void calcPulseVar() const;

        void updating() const;

    public:
        declareClassNames;
        declareObjFty(solutionWrite, controller,
                      (const flowModel &flows, const iteration &iter, const runtimeMesh &mesh,
                       const controller &cont),
                      (flows, iter, mesh, cont));

        /**
         * \brief Construct from controller and runtimeMesh.
         */
        solutionWrite(const flowModel &flows, const iteration &iter, const runtimeMesh &mesh,
                      const controller &cont);

        static uniquePtr<solutionWrite> creator(const flowModel &flows, const iteration &iter,
                                                const runtimeMesh &mesh, const controller &cont);

        /**
         * \brief Destructor.
         */
        virtual ~solutionWrite() noexcept {}

        hur_nodiscard inline const flowModel &flows() const noexcept;

        hur_nodiscard inline const iteration &iter() const noexcept;

        hur_nodiscard inline const runtimeMesh &mesh() const noexcept;

        void write() const;

        hur_nodiscard inline bool found(const string &name) const;

    private:
        void writeProfiles() const;

    protected:
        virtual void writeToFile() const = 0;

        hur_nodiscard bool samplingNow() const noexcept;

    public:
        hur_nodiscard bool writeNow() const noexcept;

        void setOutputField(const string &name, const realArray &val) const;
        void setOutputField(const string &name, realArray &&val) const;
        void setOutputField(const string &name, const string &outName, const realArray &val) const;
        void setOutputField(const string &name, const string &outName, realArray &&val) const;

        void setOutputField(const string &name, const vectorArray &val) const;
        void setOutputField(const string &name, vectorArray &&val) const;
        void setOutputField(const string &name, const string &outName,
                            const vectorArray &val) const;
        void setOutputField(const string &name, const string &outName, vectorArray &&val) const;

        void setOutputField(const string &name, const vector2DArray &val) const;
        void setOutputField(const string &name, vector2DArray &&val) const;
        void setOutputField(const string &name, const string &outName,
                            const vector2DArray &val) const;
        void setOutputField(const string &name, const string &outName, vector2DArray &&val) const;

        void setOutputField(const string &name, const tensorArray &val) const;
        void setOutputField(const string &name, tensorArray &&val) const;
        void setOutputField(const string &name, const string &outName,
                            const tensorArray &val) const;
        void setOutputField(const string &name, const string &outName, tensorArray &&val) const;

        void setOutputField(const string &name, const symmTensorArray &val) const;
        void setOutputField(const string &name, symmTensorArray &&val) const;
        void setOutputField(const string &name, const string &outName,
                            const symmTensorArray &val) const;
        void setOutputField(const string &name, const string &outName, symmTensorArray &&val) const;

    protected:
        void calcStagnationParameters() const;
        void calcViscousRatio() const;
        void calcVorticity() const;
        void calcMoleSpecies() const;

        hur_nodiscard inline const auto &rho() const noexcept;
        hur_nodiscard inline const auto &p() const noexcept;
        hur_nodiscard inline const auto &gama() const noexcept;
        hur_nodiscard inline const auto &T() const noexcept;
        hur_nodiscard inline const auto &yi() const noexcept;
        hur_nodiscard inline const auto &v() const noexcept;
    };
} // namespace OpenHurricane

#include "solutionWrite.inl"