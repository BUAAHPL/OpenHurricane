/*!
 * \file writeFieldVar.hpp
 * \brief Headers of base class of writeFieldVar.
 *        The subroutines and functions are in the <i>writeFieldVar.cpp</i> file.
 * \author Chen Zhenyi
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

#include "combustionModel.hpp"
#include "flowModel.hpp"
#include "geometryArray.hpp"
#include "iteration.hpp"
#include <set>

namespace OpenHurricane {
    class calculateFieldVar;
    /*!\brief The class of writeFieldVar.*/
    class writeFieldVar {
    protected:
        /**\brief The flow model.*/
        const flowModel &flows_;

        const iteration &iter_;

        /** \brief The name list of written variables. */
        stringList writeFieldVarList_;

        std::set<std::string> sets_;

        /** \brief If write? */
        bool isWritting_;

        /** \brief The id of this write. */
        string writeId_;

        void removeDuplicate(stdStringList &writeList) const;

        void getVarList(const controller &cont);

    private:
        mutable PtrList<cellRealArray> *yiAvePtr_;

    protected:
        inline PtrList<cellRealArray> &yiAve();

        inline bool yiSet() const;

        /**\brief Map for output variables of flow field.*/
        std::map<std::string, object *> &outFieldVarMap_;

        /**\brief Set output varaibales map. */
        void setoutVarMap(const stringList &outFieldVarList);

    public:
        declareClassNames;
        declareObjFty(writeFieldVar, controller,
                      (const flowModel &_flows, const iteration &_iter, const controller &_cont,
                       const string &_writeId, std::map<std::string, object *> &_outFieldVarMap),
                      (_flows, _iter, _cont, _writeId, _outFieldVarMap));

        // Constructors
        writeFieldVar(const flowModel &flows, const iteration &iter, const controller &cont,
                      const string &writeId, std::map<std::string, object *> &outFieldVarMap);

        static uniquePtr<writeFieldVar> creator(const flowModel &flows, const iteration &iter,
                                                const controller &cont, const string &writeId,
                                                std::map<std::string, object *> &outFieldVarMap);

        /*!\brief Destructor.*/
        virtual ~writeFieldVar() noexcept { HurDelete(yiAvePtr_); }

        virtual void updating() = 0;
        virtual void fieldOtherVarsWriting(const string &type, const combustionModel *chemtryPtr);
        virtual void fieldOtherVarsWriting(const string &type);

        void calcTotalPressure(const string &Ma, const string &pt, const realArray &p,
                               const realArray &gama, const vectorArray &v,
                               const realArray &rho) const;

        void calcTotalTemperature(const string &Ma, const string &Tt, const realArray &T,
                                  const realArray &gama, const vectorArray &v, const realArray &p,
                                  const realArray &rho) const;

        inline const runtimeMesh &mesh() const;

        inline const cellRealArray &rho() const;
        inline const cellRealArray &p() const;
        inline const cellRealArray &T() const;
        inline const cellRealArray &E() const;
        inline const cellRealArray &gama() const;
        inline const cellVectorArray &v() const;
        inline const PtrList<cellRealArray> &yi() const;

        inline void setOutputField(const string &varName, const realArray &value) const;

        inline void setOutputField(const string &varName, realArray &&value) const;

        inline void setOutputField(const string &varName, const vectorArray &value) const;

        inline void setOutputField(const string &varName, vectorArray &&value) const;

        const std::map<std::string, object *> &getoutVarMap() const;

        template <class Type, class GeometryMesh>
        Array<Type> calcTimeAveVar(const geometryArray<Type, GeometryMesh> &var,
                                   const real curretTimeGap) const {
            return var.getTimeSumPtr() / curretTimeGap;
        }

        virtual void correct() const;

        virtual void reseting() const;

        inline const string &writeId() const noexcept;

    protected:
        virtual fileName getFileName() const = 0;

        virtual bool writeNow() const = 0;

        void writeToFile() const;

    public:
        inline void write() const;
    };
} // namespace OpenHurricane

#include "writeFieldVar.inl"