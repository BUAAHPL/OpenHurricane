/*!
 * \file thermoList.cpp
 * \brief Main subroutines for thermo properties list.
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
#include "thermoList.hpp"

const OpenHurricane::real OpenHurricane::thermoList::tol_ = 1.0e-4;
const int OpenHurricane::thermoList::maxIter_ = 100;

OpenHurricane::thermoList::thermoList()
    : object("thermoList", registerTable()), species_(speciesList::nullObject()), thList_(),
      eosPtr_(nullptr), TLow_(50.0), THigh_(5000.0), globalLowT_(300.0), globalCommonT_(1000.0),
      globaHighT_(5000.0) {}

OpenHurricane::thermoList::thermoList(const speciesList &species)
    : object("thermoList", registerTable()), species_(species), thList_(), eosPtr_(nullptr),
      TLow_(50.0), THigh_(5000.0), globalLowT_(300.0), globalCommonT_(1000.0), globaHighT_(5000.0) {
}

OpenHurricane::thermoList::thermoList(const speciesList &species, const runtimeMesh &mesh)
    : object("thermoList", mesh), species_(species), thList_(), eosPtr_(nullptr), TLow_(50.0),
      THigh_(5000.0), globalLowT_(300.0), globalCommonT_(1000.0), globaHighT_(5000.0) {}

OpenHurricane::thermoList::thermoList(const speciesList &species, const controller &cont)
    : object("thermoList", registerTable()), species_(species), thList_(species.size()),
      eosPtr_(nullptr), TLow_(cont.findOrDefault("TLow", 50.0)),
      THigh_(cont.findOrDefault("THigh", 5000.0)),
      globalLowT_(cont.findOrDefault("globalLowT", 300.0)),
      globalCommonT_(cont.findOrDefault("globalCommonT", 1000.0)),
      globaHighT_(cont.findOrDefault("globalHighT", 5000.0)) {
    eosPtr_ = equationOfState::creator(species, cont.subController(equationOfState::className_));

    for (integer i = 0; i < species.size(); ++i) {
        if (cont.found(thermo::className_)) {
            if (cont.subController(thermo::className_).found(species_[i].name())) {
                thList_.set(
                    i, thermo::creator(
                           cont.subController(thermo::className_).subController(species_[i].name()),
                           *eosPtr_, i)
                           .release());
            } else {
                thList_.set(
                    i,
                    thermo::creator(cont.subController(thermo::className_), *eosPtr_, i).release());
            }
        } else {
            thList_.set(
                i, thermo::creator(cont.subController(species_[i].name()), *eosPtr_, i).release());
        }
    }
}

OpenHurricane::thermoList::thermoList(const speciesList &species, const controller &cont,
                                    const runtimeMesh &mesh)
    : object("thermoList", mesh), species_(species), thList_(species.size()), eosPtr_(nullptr),
      TLow_(cont.findOrDefault("TLow", 50.0)), THigh_(cont.findOrDefault("THigh", 5000.0)),
      globalLowT_(cont.findOrDefault("globalLowT", 300.0)),
      globalCommonT_(cont.findOrDefault("globalCommonT", 1000.0)),
      globaHighT_(cont.findOrDefault("globalHighT", 5000.0)) {
    eosPtr_ = equationOfState::creator(species, cont.subController(equationOfState::className_));

    for (integer i = 0; i < species.size(); ++i) {
        if (cont.found(thermo::className_)) {
            if (cont.subController(thermo::className_).found(species_[i].name())) {
                thList_.set(
                    i, thermo::creator(
                           cont.subController(thermo::className_).subController(species_[i].name()),
                           *eosPtr_, i)
                           .release());
            } else {
                thList_.set(
                    i,
                    thermo::creator(cont.subController(thermo::className_), *eosPtr_, i).release());
            }
        } else {
            thList_.set(
                i, thermo::creator(cont.subController(species_[i].name()), *eosPtr_, i).release());
        }
    }
}

OpenHurricane::thermoList::thermoList(const thermoList &tt)
    : object(tt), species_(tt.species_), thList_(tt.thList_), eosPtr_(nullptr), TLow_(tt.TLow_),
      THigh_(tt.THigh_), globalLowT_(tt.globalLowT_), globalCommonT_(tt.globalCommonT_),
      globaHighT_(tt.globaHighT_) {
    if (tt.eosPtr_) {
        eosPtr_ = tt.eosPtr_->clone();
    }
}

OpenHurricane::thermoList::thermoList(thermoList &&tt) noexcept
    : object(std::move(tt)), species_(tt.species_), thList_(std::move(tt.thList_)),
      eosPtr_(nullptr), TLow_(tt.TLow_), THigh_(tt.THigh_), globalLowT_(tt.globalLowT_),
      globalCommonT_(tt.globalCommonT_), globaHighT_(tt.globaHighT_) {
    if (tt.eosPtr_) {
        eosPtr_ = std::move(tt.eosPtr_);
    }
}

OpenHurricane::thermoList::~thermoList() noexcept {}

void OpenHurricane::thermoList::setEos(const controller &cont) {
    eosPtr_.clear();

    eosPtr_ = equationOfState::creator(species_, cont);
    for (integer i = 0; i < thList_.size(); ++i) {
        thList_[i].setEos(eosPtr_.get());
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::limit(const real T) const {
    if (T < TLow_ || T > THigh_) {
#ifdef HUR_DEBUG
        std::string errMsg = " Exceed temperature range of thermo table ";
        errMsg += toString(TLow_);
        errMsg += " ~ ";
        errMsg += toString(THigh_);
        errMsg += ";  T = ";
        errMsg += toString(T);
        checkWarningStr(errMsg);
#endif // HUR_DEBUG
        return min(max(T, TLow_), THigh_);
    } else {
        return T;
    }
}

OpenHurricane::real OpenHurricane::thermoList::limit(const real T, integer &iFlag) const {
    if (T < TLow_) {
        iFlag = 0;
        return TLow_;
    } else if (T > THigh_) {
        iFlag = 2;
        return THigh_;
    }

    iFlag = 1;
    return T;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::cp0(const real T, const realArray &yi) const {
#ifdef HUR_DEBUG

    if (yi.size() != species_.size()) {
        std::string errMsg;
        errMsg = "The size of species list: ";
        errMsg += toString(species_.size());
        errMsg += " is not equal to Yi(";
        errMsg += toString(yi.size());
        errMsg += ")";
        errorAbortStr(errMsg);
    }
#endif // HUR_DEBUG
    if (species_.size() == 1) {
        return thList_[0].cp0(T);
    }
    real cpm0 = 0.0;
    for (integer i = 0; i < yi.size(); ++i) {
        /* Pout.setReal();
         Pout << species_.W(i) << std::endl;*/
        cpm0 += thList_[i].cp0(T) * yi[i];
        /* Pout.unsetReal();*/
    }

    return cpm0;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::cp0(const real T,
                                                          const PtrList<cellRealArray> &yi,
                                                          const integer cellI) const {
    if (species_.size() == 1) {
        return thList_[0].cp0(T);
    }
    real cpm0 = 0.0;
    for (integer i = 0; i < yi.size(); ++i) {
        cpm0 += thList_[i].cp0(T) * yi[i][cellI];
    }
    return cpm0;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::Cp0(const real T, const realArray &xi) const {
#ifdef HUR_DEBUG

    if (xi.size() != species_.size()) {
        std::string errMsg;
        errMsg = "The size of species list: ";
        errMsg += toString(species_.size());
        errMsg += " is not equal to Xi(";
        errMsg += toString(xi.size());
        errMsg += ")";
        errorAbortStr(errMsg);
    }
#endif // HUR_DEBUG
    if (species_.size() == 1) {
        return thList_[0].Cp0(T);
    }
    real Cpm0 = 0.0;
    for (integer i = 0; i < xi.size(); ++i) {
        Cpm0 += thList_[i].Cp0(T) * xi[i];
    }
    return Cpm0;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::cv0(const real T, const realArray &yi) const {
#ifdef HUR_DEBUG

    if (yi.size() != species_.size()) {
        std::string errMsg;
        errMsg = "The size of species list: ";
        errMsg += toString(species_.size());
        errMsg += " is not equal to Yi(";
        errMsg += toString(yi.size());
        errMsg += ")";
        errorAbortStr(errMsg);
    }

#endif // HUR_DEBUG
    if (species_.size() == 1) {
        return thList_[0].cv0(T);
    }
    real cvm0 = 0.0;
    for (integer i = 0; i < yi.size(); ++i) {
        cvm0 += thList_[i].cv0(T) * yi[i];
    }
    return cvm0;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::cv0(const real T,
                                                          const PtrList<cellRealArray> &yi,
                                                          const integer cellI) const {
    if (species_.size() == 1) {
        return thList_[0].cv0(T);
    }
    real cvm0 = 0.0;
    for (integer i = 0; i < yi.size(); ++i) {
        cvm0 += thList_[i].cv0(T) * yi[i][cellI];
    }
    return cvm0;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::ha0(const real T, const realArray &yi) const {
#ifdef HUR_DEBUG

    if (yi.size() != species_.size()) {
        std::string errMsg;
        errMsg = "The size of species list: ";
        errMsg += toString(species_.size());
        errMsg += " is not equal to Yi(";
        errMsg += toString(yi.size());
        errMsg += ")";
        errorAbortStr(errMsg);
    }
#endif // HUR_DEBUG
    if (species_.size() == 1) {
        return thList_[0].ha0(T);
    }
    real ham0 = 0.0;
    for (integer i = 0; i < yi.size(); ++i) {
        ham0 += thList_[i].ha0(T) * yi[i];
    }
    return ham0;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::ha0(const real T,
                                                          const PtrList<cellRealArray> &yi,
                                                          const integer cellI) const {
    if (species_.size() == 1) {
        return thList_[0].ha0(T);
    }
    real ham0 = 0.0;
    for (integer i = 0; i < yi.size(); ++i) {
        ham0 += thList_[i].ha0(T) * yi[i][cellI];
    }
    return ham0;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::thermoList::ha0(const real T) const {
    realArray tempHa0(species_.size(), Zero);

    for (integer i = 0; i < tempHa0.size(); ++i) {
        tempHa0[i] = thList_[i].ha0(T);
    }

    return tempHa0;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::thermoList::hc() const {
    realArray tempHc0(species_.size(), Zero);

    for (integer i = 0; i < tempHc0.size(); ++i) {
        tempHc0[i] = thList_[i].hc();
    }

    return tempHc0;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::hc(const realArray &yi) const {
#ifdef HUR_DEBUG

    if (yi.size() != species_.size()) {
        std::string errMsg;
        errMsg = "The size of species list: ";
        errMsg += toString(species_.size());
        errMsg += " is not equal to Yi(";
        errMsg += toString(yi.size());
        errMsg += ")";
        errorAbortStr(errMsg);
    }

#endif // HUR_DEBUG
    if (species_.size() == 1) {
        return thList_[0].hc();
    }
    real hc0 = 0.0;
    for (integer i = 0; i < yi.size(); ++i) {
        hc0 += thList_[i].hc() * yi[i];
    }
    return hc0;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::hc(const PtrList<cellRealArray> &yi,
                                                         const integer cellI) const {
    if (species_.size() == 1) {
        return thList_[0].hc();
    }
    real hc0 = 0.0;
    for (integer i = 0; i < yi.size(); ++i) {
        hc0 += thList_[i].hc() * yi[i][cellI];
    }
    return hc0;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::s_p(const real pm, const real T,
                                                          const realArray &yi) const {
#ifdef HUR_DEBUG

    if (yi.size() != species_.size()) {
        std::string errMsg;
        errMsg = "The size of species list: ";
        errMsg += toString(species_.size());
        errMsg += " is not equal to Yi(";
        errMsg += toString(yi.size());
        errMsg += ")";
        errorAbortStr(errMsg);
    }

#endif // HUR_DEBUG
    if (species_.size() == 1) {
        return thList_[0].s_p(pm, T);
    }
    realArray xi = species_.Yi2Xi(yi);
    real sm = Zero;
    for (integer i = 0; i < yi.size(); ++i) {
        sm += thList_[i].s_p(xi[i] * pm, T) * yi[i];
    }
    return sm;
}

hur_nodiscard OpenHurricane::real OpenHurricane::thermoList::s_p(const real pm, const real T,
                                                          const PtrList<cellRealArray> &yi,
                                                          const integer cellI) const {
    if (species_.size() == 1) {
        return thList_[0].s_p(pm, T);
    }
    realArray xi = species_.Yi2Xi(yi, cellI);
    real sm = Zero;
    for (integer i = 0; i < xi.size(); ++i) {
        sm += thList_[i].s_p(xi[i] * pm, T) * yi[i][cellI];
    }
    return sm;
}

OpenHurricane::real OpenHurricane::thermoList::inteCp0dT(const real T1, const real T2,
                                                  const realArray &yi) const {
    if (species_.size() == 1) {
        return thList_[0].inteCp0dT(T1, T2);
    }

    real sm = Zero;
    for (integer i = 0; i < species_.size(); ++i) {
        sm += thList_[i].inteCp0dT(T1, T2) * yi[i];
    }
    return sm;
}

OpenHurricane::real OpenHurricane::thermoList::inteCp0dT(const real T1, const real T2,
                                                  const PtrList<cellRealArray> &yi,
                                                  const integer cellI) const {
    if (species_.size() == 1) {
        return thList_[0].inteCp0dT(T1, T2);
    }
    real sm = Zero;
    for (integer i = 0; i < species_.size(); ++i) {
        sm += thList_[i].inteCp0dT(T1, T2) * yi[i][cellI];
    }
    return sm;
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::thermoList::g0(const real T) const {
    realArray tempg0(species_.size(), Zero);

    for (integer i = 0; i < tempg0.size(); ++i) {
        tempg0[i] = thList_[i].g0(T);
    }

    return tempg0;
}

void OpenHurricane::thermoList::g0(const real T, realArray &gf) const {
    if (gf.size() != species_.size()) {
        gf.resize(species_.size(), Zero);
    }
    for (integer i = 0; i < gf.size(); ++i) {
        gf[i] = thList_[i].g0(T);
    }
}

hur_nodiscard OpenHurricane::realArray OpenHurricane::thermoList::G0(const real T) const {
    realArray tempg0(species_.size(), Zero);

    for (integer i = 0; i < tempg0.size(); ++i) {
        tempg0[i] = thList_[i].G0(T);
    }

    return tempg0;
}

void OpenHurricane::thermoList::G0(const real T, realArray &Gf) const {
    if (Gf.size() != species_.size()) {
        Gf.resize(species_.size(), Zero);
    }
    for (integer i = 0; i < Gf.size(); ++i) {
        Gf[i] = thList_[i].G0(T);
    }
}

void OpenHurricane::thermoList::DG0DT(const real T, realArray &DG0DTf) const {
    if (DG0DTf.size() != species_.size()) {
        DG0DTf.resize(species_.size(), Zero);
    }
    for (integer i = 0; i < DG0DTf.size(); ++i) {
        DG0DTf[i] = thList_[i].DG0DT(T);
    }
}

OpenHurricane::real OpenHurricane::thermoList::T(
    real f, real p, real T0, integer &iFlag, const PtrList<cellRealArray> &yi, const integer cellI,
    real (thermoList::*F)(const real, const real, const PtrList<cellRealArray> &, const integer)
        const,
    real (thermoList::*dFdT)(const real, const real, const PtrList<cellRealArray> &, const integer)
        const,
    real (thermoList::*limit)(const real, integer &) const) const {
    real Test = T0;
    real Tnew = T0;
    real Ttol = T0 * tol_;
    int iter = 0;

    do {
        Test = Tnew;
        Tnew = (this->*limit)(
            Test - ((this->*F)(p, Test, yi, cellI) - f) / (this->*dFdT)(p, Test, yi, cellI), iFlag);

        if (iter++ > maxIter_) {
#ifdef HUR_DEBUG
            PLWarning("Maximum number of iterations exceeded: %d", maxIter_);
#else
            if (report) {
                LWarning("Maximum number of iterations exceeded: %d", maxIter_);
            }
#endif // HUR_DEBUG
            break;
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}

OpenHurricane::real
OpenHurricane::thermoList::T(real f, real p, real T0, integer &iFlag, const realArray &yi,
                          real (thermoList::*F)(const real, const real, const realArray &) const,
                          real (thermoList::*dFdT)(const real, const real, const realArray &)
                              const,
                          real (thermoList::*limit)(const real, integer &) const) const {
    real Test = T0;
    real Tnew = T0;
    real Ttol = T0 * tol_;
    int iter = 0;

    do {
        Test = Tnew;
        Tnew = (this->*limit)(Test - ((this->*F)(p, Test, yi) - f) / (this->*dFdT)(p, Test, yi),
                              iFlag);

        if (iter++ > maxIter_) {
#ifdef HUR_DEBUG
            PLWarning("Maximum number of iterations exceeded: %d", maxIter_);
#else
            if (report) {
                LWarning("Maximum number of iterations exceeded: %d", maxIter_);
            }
#endif // HUR_DEBUG           
            break;
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}
