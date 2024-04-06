/*!
 * \file chemistrySource.cpp
 * \brief Main subroutines for chemistry source.
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
#include "chemistrySource.hpp"
#include <functional>

 namespace OpenHurricane{
	createClassNameStr(chemistrySource,"chemistrySource");
}

namespace OpenHurricane {
     createObjFty(chemistrySource,controller);

} // namespace OpenHurricane

OpenHurricane::chemistrySource::chemistrySource(flowModel &flows, const controller &cont)
    : flows_(flows), mesh_(flows.mesh()), yi_(flows.mixtures().Yi()),
      reactions_(flows.mixtures().reactions()), species_(flows.mixtures().species()),
      nsp_(flows.mixtures().species().size()), nSpc_(flows.mixtures().species().size()),
      nrc_(flows.mixtures().reactions().size()), nRact_(flows.mixtures().reactions().size()),
      isReacting_(!flows.mixtures().noReaction()), simplifiedJacobians_(false), ci_(), yiTmp_(),
      dcdt_(), RiPtr_(nullptr),
#ifdef TEST_PROCESS_TIME
      odeCountIter_(0),
#endif
      dtInitPtr_(nullptr), dtMax_(cont.findOrDefault<real>("chemicalMaxStep", large)), hea0_(0),
      T_(0), p0_(0), p_(0), rho0_(0), rho_(0),
      reactionFlowType_(ReactionFlowType::ConstantPressure),
#ifdef TEST_PROCESS_TIME
      calcTime_(0), reducionTime_(0),
#endif
      nSPInEveryCellPtr_(nullptr), nRCInEveryCellPtr_(nullptr), diagJacPerReacPtr_(nullptr),
      TemperatureAlpha_(0.2), TemperatureBeta_(0.25), TemperatureThreshold_(200),
      cellSourceCalTimePtr_(nullptr) {
    if (cont.found("ReactionFlowType")) {
        string rw = cont.findWord("ReactionFlowType");
        if (rw == "ConstantPressure") {
            reactionFlowType_ = ReactionFlowType::ConstantPressure;
        } else if (rw == "ConstantVolume") {
            reactionFlowType_ = ReactionFlowType::ConstantVolume;
        } else if (rw == "IsolatedSystem") {
            reactionFlowType_ = ReactionFlowType::IsolatedSystem;
        } else {
            errorAbortStr(("Unknown reaction flow type: " + rw + " in " + cont.name()));
        }
    }
    if (isReacting_) {
        ci_.resize(nsp_, Zero);
        yiTmp_.resize(nsp_, Zero);
        dcdt_.resize(nsp_, Zero);
    }
    simplifiedJacobians_ = controllerSwitch(cont)("simplifiedJacobians", simplifiedJacobians_);
    if (cont.found("solveOptions")) {
        const auto &solOptCont = cont.subController("solveOptions");
        TemperatureThreshold_ =
            solOptCont.findOrDefault<real>("TemperatureThreshold", TemperatureThreshold_);
        TemperatureAlpha_ = solOptCont.findOrDefault<real>("TemperatureAlpha", TemperatureAlpha_);
        TemperatureBeta_ = solOptCont.findOrDefault<real>("TemperatureBeta", TemperatureBeta_);
    }
}

OpenHurricane::uniquePtr<OpenHurricane::chemistrySource>
OpenHurricane::chemistrySource::creator(flowModel &flows, const controller &cont) {
    string solverType = cont.findWord(chemistrySource::className_);

    
	defineInObjCreator(chemistrySource,static_cast<std::string>(solverType),controller,(flows, cont));
}

OpenHurricane::chemistrySource::~chemistrySource() noexcept {
    HurDelete(RiPtr_);
    HurDelete(diagJacPerReacPtr_);
    HurDelete(dtInitPtr_);
    HurDelete(nSPInEveryCellPtr_);
    HurDelete(nRCInEveryCellPtr_);

#ifdef TEST_PROCESS_TIME
    Pout.setReal();
    Pout << "    The time consumed by calculation is " << calcTime_ << " [s]" << std::endl;
    Pout.unsetReal();
#endif

    HurDelete(cellSourceCalTimePtr_);
}

//OpenHurricane::real OpenHurricane::chemistrySource::qj(const real p, const real T, const realArray& c, const integer j)
//{
//	real qf = 0.0;
//	real qr = 0.0;
//	return qfrj(p, T, c, j, qf, qr);
//}

OpenHurricane::real OpenHurricane::chemistrySource::qj(const real p, const real T, const realArray &c,
                                               const integer j, realArray &diagJac,
                                               const bool calcLastSpeciesJac) {
    real qf = 0.0;
    real qr = 0.0;
    return qfrj(p, T, c, j, qf, qr, diagJac, calcLastSpeciesJac);
}

//OpenHurricane::real OpenHurricane::chemistrySource::qfrj
//(
//	const real p,
//	const real T,
//	const realArray& c,
//	const integer j,
//	real& qfj,
//	real& qrj
//)
//{
//	Reactions& reaction = reactions()[j];
//	real kf = reaction.kf(p, T, c);
//	real kr = reaction.kr(kf, p, T, c);
//
//	const auto& fc = reaction.forwardCoeffs();
//	real qf = kf;
//	for (integer i = 0; i < fc.size(); ++i)
//	{
//		const integer spi = fc[i].index_;
//		qf *= pow(c[spi], fc[i].order_);
//	}
//
//	const auto& bc = reaction.backwardCoeffs();
//	real qr = kr;
//	for (integer i = 0; i < bc.size(); ++i)
//	{
//		const integer spi = bc[i].index_;
//		qr *= pow(c[spi], bc[i].order_);
//	}
//	qfj = qf;
//	qrj = qr;
//	return qf - qr;
//}

OpenHurricane::real OpenHurricane::chemistrySource::qfrj(const real p, const real T, const realArray &c,
                                                 const integer j, real &qfj, real &qrj,
                                                 realArray &diagJac,
                                                 const bool calcLastSpeciesJac) {
    auto &reaction = reactions()[j];
    real kf = reaction.kf(p, T, c);
    real kr = reaction.kr(kf, p, T, c);

    diagJac = Zero;
    const auto &fc = reaction.forwardCoeffs();
    real qf = kf;
    for (integer i = 0; i < fc.size(); ++i) {
        const integer spi = fc[i].index_;
        qf *= pow(c[spi], fc[i].order_);
    }

    for (integer i = 0; i < fc.size(); ++i) {
        const integer spi = fc[i].index_;
        real pf = kf;
        for (integer k = 0; k < fc.size(); ++k) {
            const integer sk = fc[k].index_;
            const real orderk = fc[k].order_;

            if (i == k) {
                if (orderk < 1.0) {
                    if (c[sk] > tiny) {
                        pf *= orderk * pow(c[sk] + veryTiny, orderk - real(1.0));
                    } else {
                        pf = real(0.0);
                    }
                } else if (orderk == 1.0) {
                    continue;
                } else {
                    pf *= orderk * pow(c[sk], orderk - real(1.0));
                }
            } else {
                pf *= pow(c[sk], orderk);
            }
        }
        if (calcLastSpeciesJac) {
            diagJac[spi] += pf;
        } else {
            if (spi != nsp() - 1) {
                diagJac[spi] += pf;
            } else //spi=Ns
            {
                for (integer jj = 0; jj < nsp() - 1; ++jj) {
                    diagJac[jj] -= species().W(jj) / species().W(spi) * pf;
                }
            }
        }
    }

    const auto &bc = reaction.backwardCoeffs();
    real qr = kr;
    for (integer i = 0; i < bc.size(); ++i) {
        const integer spi = bc[i].index_;
        qr *= pow(c[spi], bc[i].order_);
    }

    for (integer i = 0; i < bc.size(); ++i) {
        const integer spi = bc[i].index_;
        real pr = kr;
        for (integer k = 0; k < bc.size(); ++k) {
            const integer sk = bc[k].index_;
            const real orderk = bc[k].order_;

            if (i == k) {
                if (orderk < 1.0) {
                    if (c[sk] > tiny) {
                        pr *= orderk * pow(c[sk] + veryTiny, orderk - real(1.0));
                    } else {
                        pr = real(0.0);
                    }
                } else if (orderk == 1.0) {
                    continue;
                } else {
                    pr *= orderk * pow(c[sk], orderk - real(1.0));
                }
            } else {
                pr *= pow(c[sk], orderk);
            }
        }
        if (calcLastSpeciesJac) {
            diagJac[spi] -= pr;
        } else {
            if (spi != nsp() - 1) {
                diagJac[spi] -= pr;
            } else //spi=Ns
            {
                for (integer jj = 0; jj < nsp() - 1; ++jj) {
                    diagJac[jj] += species().W(jj) / species().W(spi) * pr;
                }
            }
        }
    }

    qfj = qf;
    qrj = qr;

    return qf - qr;
}

OpenHurricane::realArray OpenHurricane::chemistrySource::omega(const real p, const real T,
                                                       const realArray &c) {
    reactions().updateG0(T);
    realArray dcdt(nEqns(), Zero);
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        reactions()[ri].q(p, T, c, dcdt);
        /*real qri = qj(p, T, c, ri);
        for (integer i = 0; i < reactions()[ri].forwardCoeffs().size(); ++i)
        {
                integer index = reactions()[ri].forwardCoeffs()[i].index_;
                real stoic = reactions()[ri].forwardCoeffs()[i].stoichCoeff_;
                if (index < nEqns())
                {
                        dcdt[index] -= stoic * qri;
                }
        }
        for (integer i = 0; i < reactions()[ri].backwardCoeffs().size(); ++i)
        {
                integer index = reactions()[ri].backwardCoeffs()[i].index_;
                real stoic = reactions()[ri].backwardCoeffs()[i].stoichCoeff_;
                if (index < nEqns())
                {
                        dcdt[index] += stoic * qri;
                }
        }*/
    }

    return dcdt;
}

void OpenHurricane::chemistrySource::omega(const real p, const real T, const realArray &c,
                                       realArray &w) {
    reactions().updateG0(T);
    if (w.size() < nsp()) {
        w.resize(nsp());
    }
    w = Zero;
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        reactions()[ri].q(p, T, c, w);
        /*real qri = qj(p, T, c, ri);
        for (integer i = 0; i < reactions()[ri].forwardCoeffs().size(); ++i)
        {
                integer index = reactions()[ri].forwardCoeffs()[i].index_;
                real stoic = reactions()[ri].forwardCoeffs()[i].stoichCoeff_;
                w[index] -= stoic * qri;

        }
        for (integer i = 0; i < reactions()[ri].backwardCoeffs().size(); ++i)
        {
                integer index = reactions()[ri].backwardCoeffs()[i].index_;
                real stoic = reactions()[ri].backwardCoeffs()[i].stoichCoeff_;
                w[index] += stoic * qri;
        }*/
    }
}

void OpenHurricane::chemistrySource::omegaFullJacobian(const real p, const real T, const realArray &c,
                                                   realArray &w, realSquareMatrix &dwdy) {
    reactions().updateG0(T);
    w = Zero;
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        const auto &rr = reactions()[ri];
        real qfr, qf, qr, kf, kr;
        qfr = rr.q(p, T, c, w, qf, qr, kf, kr);
        rr.DqDci(kf, kr, qfr, p, T, c, dwdy);
        rr.DqDT(kf, kr, qfr, p, T, c, dwdy);
    }
}

void OpenHurricane::chemistrySource::omega(const real rho, const real p, const real T,
                                       const realArray &c, const realArray &yi, realArray &w,
                                       realArray &diagJac) {
    reactions().updateG0(T);
    if (w.size() < nsp()) {
        w.resize(nsp());
    }
    if (diagJac.size() < nsp()) {
        diagJac.resize(nsp());
    }
    w = Zero;
    diagJac = Zero;
    auto &diagSource = diagJacPerReac();
    if (diagSource.size() != nsp()) {
        diagSource.resize(nsp());
    }
    if (simplifiedJacobians_) {
        for (integer ri = 0; ri < reactions().size(); ++ri) {
            diagSource = Zero;
            real qri = qj(p, T, c, ri, diagSource);

            const auto &fc = reactions()[ri].forwardCoeffs();
            for (integer i = 0; i < fc.size(); ++i) {
                const integer spi = fc[i].index_;
                const real stoich = fc[i].stoichCoeff_;
                /*if (spi < nsp() - 1)
                {*/
                w[spi] -= stoich * qri;
                diagJac[spi] -= stoich * diagSource[spi];
                //}
            }

            // vr
            const auto &bc = reactions()[ri].backwardCoeffs();
            for (integer i = 0; i < bc.size(); ++i) {
                const integer spi = bc[i].index_;
                const real stoich = bc[i].stoichCoeff_;
                /*if (spi < nsp() - 1)
                {*/
                w[spi] += stoich * qri;
                diagJac[spi] += stoich * diagSource[spi];
                /*}*/
            }
        }
    } else {
        diagSource = Zero;
        for (integer ri = 0; ri < reactions().size(); ++ri) {
            real kf, kr;
            real qfr = reactions()[ri].q(p, T, c, w, diagJac, true, kf, kr);
            reactions()[ri].DqDT(kf, kr, qfr, p, T, c, diagSource);
        }
        const real cv = therm().cv_p(p, T, yi);
        for (integer i = 0; i < nsp_; ++i) {
            diagJac[i] +=
                diagSource[i] * therm()[i].Wi() * (-therm()[i].ea(rho, p, T)) / (rho * cv);
        }
    }
}

void OpenHurricane::chemistrySource::omegaCoupled(const real p, const real T, const realArray &c,
                                              realArray &w) {
    reactions().updateG0(T);
    if (w.size() < nsp()) {
        w.resize(nsp());
    }
    w = Zero;
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        const auto &reac = reactions()[ri];
        real qf, qr;
        real qri = reac.q(p, T, c, qf, qr);
        //real qri = qj(p, T, c, ri);
        for (integer i = 0; i < reactions()[ri].forwardCoeffs().size(); ++i) {
            integer index = reactions()[ri].forwardCoeffs()[i].index_;
            real stoic = reactions()[ri].forwardCoeffs()[i].stoichCoeff_;
            /*if (index < nsp() - 1)
            {*/
            w[index] -= stoic * qri;
            //}
        }
        for (integer i = 0; i < reactions()[ri].backwardCoeffs().size(); ++i) {
            integer index = reactions()[ri].backwardCoeffs()[i].index_;
            real stoic = reactions()[ri].backwardCoeffs()[i].stoichCoeff_;
            /*if (index < nsp() - 1)
            {*/
            w[index] += stoic * qri;
            /*}*/
        }
    }
}

void OpenHurricane::chemistrySource::omegaCoupled(const real rho, const real p, const real T,
                                              const realArray &c, const realArray &yi, realArray &w,
                                              realArray &diagJac) {
    reactions().updateG0(T);
    if (w.size() < nsp()) {
        w.resize(nsp());
    }
    if (diagJac.size() < nsp() - 1) {
        diagJac.resize(nsp() - 1);
    }
    w = Zero;
    diagJac = Zero;
    auto &diagSource = diagJacPerReac();
    if (diagSource.size() != nsp()) {
        diagSource.resize(nsp());
    }
    if (simplifiedJacobians_) {
        for (integer ri = 0; ri < reactions().size(); ++ri) {
            diagSource = Zero;
            real qri = qj(p, T, c, ri, diagSource);

            const auto &fc = reactions()[ri].forwardCoeffs();
            for (integer i = 0; i < fc.size(); ++i) {
                const integer spi = fc[i].index_;
                const real stoich = fc[i].stoichCoeff_;
                w[spi] -= stoich * qri;
                if (spi < nsp() - 1) {
                    diagJac[spi] -= stoich * diagSource[spi];
                }
            }

            // vr
            const auto &bc = reactions()[ri].backwardCoeffs();
            for (integer i = 0; i < bc.size(); ++i) {
                const integer spi = bc[i].index_;
                const real stoich = bc[i].stoichCoeff_;
                w[spi] += stoich * qri;
                if (spi < nsp() - 1) {
                    diagJac[spi] += stoich * diagSource[spi];
                }
            }
        }
    } else {
        diagSource = Zero;
        for (integer ri = 0; ri < reactions().size(); ++ri) {
            real kf, kr;
            real qfr = reactions()[ri].q(p, T, c, w, diagJac, false, kf, kr);
            reactions()[ri].DqDT(kf, kr, qfr, p, T, c, diagSource);
        }
        const real cv = therm().cv_p(p, T, yi);
        const real ens = therm()[nsp_ - 1].ea(rho, p, T);
        const real betae = inv(rho * cv);
        for (integer i = 0; i < nsp_ - 1; ++i) {
            //diagJac[i] += diagSource[i] * therm()[i].Wi() * (ens - therm()[i].ea(rho, p, T)) / (rho * cv);
            diagJac[i] +=
                diagSource[i] * therm()[i].Wi() * (ens - therm()[i].ea(rho, p, T)) * betae;
        }
    }
}

OpenHurricane::real OpenHurricane::chemistrySource::omegaCoupledAndtc(const real p, const real T,
                                                              const realArray &c,
                                                              realArray &omega) {
    reactions().updateG0(T);
    if (omega.size() < nsp()) {
        omega.resize(nsp());
    }
    omega = Zero;
    //real lastSpcOmega = Zero;
    real ciSum = Zero;

    for (integer i = 0; i < nsp(); ++i) {
        ciSum += c[i];
    }
    real ts = Zero;
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        const auto &reac = reactions()[ri];
        real qf, qr;
        real qri = reac.q(p, T, c, qf, qr);
        //real qri = qj(p, T, c, ri);
        for (integer i = 0; i < reactions()[ri].forwardCoeffs().size(); ++i) {
            integer index = reactions()[ri].forwardCoeffs()[i].index_;
            real stoic = reactions()[ri].forwardCoeffs()[i].stoichCoeff_;
            /*if (index < nsp() - 1)
            {*/
            omega[index] -= stoic * qri;
            //}
            /*else
            {
                    lastSpcOmega -= stoic * qri;
            }*/
        }
        for (integer i = 0; i < reactions()[ri].backwardCoeffs().size(); ++i) {
            integer index = reactions()[ri].backwardCoeffs()[i].index_;
            real stoic = reactions()[ri].backwardCoeffs()[i].stoichCoeff_;

            /*if (index < nsp() - 1)
            {*/
            omega[index] += stoic * qri;
            //}
            /*else
            {
                    lastSpcOmega += stoic * qri;
            }*/
            ts += stoic * qf;
        }
    }
    /*realArray tci(nsp(), Zero);

    tci[nsp() - 1] = c[nsp() - 1] / mag(lastSpcOmega);

    if (isnan(tci[nsp() - 1]) || isinf(tci[nsp() - 1]))
    {
            tci[nsp() - 1] = 0.0;
    }
    for (integer i = 0; i < nsp() - 1; ++i)
    {
            tci[i] = c[i] / mag(omega[i]);
            if (isnan(tci[i]) || isinf(tci[i]))
            {
                    tci[i] = 0.0;
            }
    }
    real tmp = tci[0];
    for (integer i = 1; i < nsp(); ++i)
    {
            tmp = max(tmp, tci[i]);
    }
    return tmp;*/

    /*real tcc = c[nsp() - 1] / mag(lastSpcOmega);
    if (isnan(tcc) || isinf(tcc))
    {
            tcc = 0;
    }
    for (integer i = 0; i < nsp() - 1; ++i)
    {
            real tcii = c[i] / mag(omega[i]);
            if (isnan(tcii) || isinf(tcii))
            {
                    tcii = 0;
            }
            tcc = max(tcc, tcii);
    }
    return tcc;*/
    ts = real(reactions().size()) * ciSum / ts;
    if (isnan(ts) || isinf(ts)) {
        return large;
    }
    return min(ts, large);
}

OpenHurricane::real OpenHurricane::chemistrySource::omegaCoupledAndtc(const real rho, const real p,
                                                              const real T, const realArray &c,
                                                              const realArray &yi, realArray &omega,
                                                              realArray &diagJac) {
    reactions().updateG0(T);
    if (omega.size() < nsp()) {
        omega.resize(nsp());
    }
    if (diagJac.size() < nsp() - 1) {
        diagJac.resize(nsp() - 1);
    }
    omega = Zero;
    diagJac = Zero;
    auto &diagSource = diagJacPerReac();
    if (diagSource.size() < nsp()) {
        diagSource.resize(nsp());
    }
    //real lastSpcOmega = Zero;
    real ciSum = Zero;

    for (integer i = 0; i < nsp(); ++i) {
        ciSum += c[i];
    }
    real ts = Zero;
    diagSource = Zero;
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        real qf = 0;
        real qr = 0;
        real kf, kr;
        real qfr = reactions()[ri].q(p, T, c, omega, qf, qr, kf, kr);
        reactions()[ri].DqDci(kf, kr, qfr, p, T, c, diagJac, false);
        reactions()[ri].DqDT(kf, kr, qfr, p, T, c, diagSource);

        // vr
        const auto &bc = reactions()[ri].backwardCoeffs();
        for (integer i = 0; i < bc.size(); ++i) {
            const integer spi = bc[i].index_;
            const real stoich = bc[i].stoichCoeff_;
            /*if (spi < nsp() - 1)
            {
                    omega[spi] += stoich * qri;
                    diagJac[spi] += stoich * diagSource[spi];
            }*/
            /*else
            {
                    lastSpcOmega += stoich * qri;
            }*/
            ts += stoich * qf;
        }
    }
    /*real tcc = c[nsp() - 1] / mag(lastSpcOmega);
    if (isnan(tcc) || isinf(tcc))
    {
            tcc = 0;
    }
    for (integer i = 0; i < nsp() - 1; ++i)
    {
            real tcii = c[i] / mag(omega[i]);
            if (isnan(tcii) || isinf(tcii))
            {
                    tcii = 0;
            }
            tcc = max(tcc, tcii);
    }
    return tcc;*/
    const real cv = therm().cv_p(p, T, yi);
    const real ens = therm()[nsp_ - 1].ea(rho, p, T);
    const real betae = inv(rho * cv);
    for (integer i = 0; i < nsp_ - 1; ++i) {
        //diagJac += diagSource[i] * therm()[i].Wi() * (ens - therm()[i].ea(rho, p, T)) / (rho * cv);
        diagJac[i] += diagSource[i] * therm()[i].Wi() * (ens - therm()[i].ea(rho, p, T)) * betae;
    }
    ts = real(reactions().size()) * ciSum / ts;
    if (isnan(ts) || isinf(ts)) {
        return large;
    }
    return min(ts, large);
}

void OpenHurricane::chemistrySource::timescale(const real p, const real T, const realArray &c,
                                           realArray &ts) {
    reactions().updateG0(T);
    realArray diagSource(species_.size(), Zero);
    ts = Zero;
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        diagSource = Zero;
        real qri = qj(p, T, c, ri, diagSource);

        const auto &fc = reactions()[ri].forwardCoeffs();
        for (integer i = 0; i < fc.size(); ++i) {
            const integer spi = fc[i].index_;
            const real stoich = fc[i].stoichCoeff_;

            ts[spi] -= stoich * diagSource[spi];
        }

        // vr
        const auto &bc = reactions()[ri].backwardCoeffs();
        for (integer i = 0; i < bc.size(); ++i) {
            const integer spi = bc[i].index_;
            const real stoich = bc[i].stoichCoeff_;

            ts[spi] += stoich * diagSource[spi];
        }
    }
    for (integer i = 0; i < ts.size(); ++i) {
        ts[i] = mag(ts[i]);
        if (ts[i] <= tiny) {
            ts[i] = large;
        } else {
            ts[i] = inv(ts[i]);
        }

        ts[i] = min(ts[i], large);
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::chemistrySource::tcSFR(const real p, const real T,
                                                                const realArray &c) {
    realArray ts(species_.size(), Zero);
    timescale(p, T, c, ts);

    real maxTc = Zero;
    for (integer i = 0; i < ts.size(); ++i) {
        if (ts[i] >= 1e2) {
            continue;
        }
        maxTc = max(maxTc, ts[i]);
    }

    if (maxTc == 0) {
        maxTc = 1e2;
    }
    return maxTc;
}

void OpenHurricane::chemistrySource::timescale2(const real rho, const real p, const real T,
                                            const realArray &c, const realArray &yyi,
                                            realArray &ts) {
    realArray dqdc(nsp(), Zero);
    reactions().updateG0(T);
    ts = Zero;
    omega(rho, p, T, c, yyi, dcdt_, dqdc);

    for (integer i = 0; i < ts.size(); ++i) {
        ts[i] = 1.0 / max(mag(dqdc[i]), tiny);
        ts[i] = min(ts[i], large);
    }
}

void OpenHurricane::chemistrySource::timescale2(const real p, const real T, const realArray &c,
                                            realArray &ts) {
    realArray dqdc(nsp(), Zero);
    reactions().updateG0(T);
    ts = Zero;

    for (integer i = 0; i < this->nsp_; ++i) {
        this->ci_[i] = c[i];
    }
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        real qf = 0.0;
        real qr = 0.0;
        dqdc = Zero;
        real kf, kr;
        real qfr = reactions()[ri].q(p, T, this->ci_, qf, qr, kf, kr);
        reactions()[ri].DqDci(kf, kr, qfr, p, T, this->ci_, ts, true);
    }
    for (integer i = 0; i < ts.size(); ++i) {
        if (ts[i] >= 0.0) {
            ts[i] = large;
        } else {
            ts[i] = -1.0 / ts[i];
        }
    }
}

void OpenHurricane::chemistrySource::timescaleReacMode(const real p, const real T, const realArray &c,
                                                   realArray &ts) {
    realArray cT(nsp() + 1);
    realArray dqdc(nsp() + 1);
    realSquareMatrix Jac(nsp() + 1);
    reactions().updateG0(T);
    ts = Zero;
    for (integer i = 0; i < this->nsp_; ++i) {
        cT[i] = c[i];
    }
    cT[nsp()] = T;
    p_ = p;
    jacobian(0, cT, dqdc, Jac);

    auto eige = Jac.eigenvalues();

    for (integer i = 0; i < ts.size(); ++i) {
        auto lambdam = eige[i].real();
        if (mag(lambdam) < veryTiny) {
            ts[i] = large;
        } else {
            ts[i] = 1.0 / mag(lambdam);
        }

        ts[i] = min(real(200), ts[i]);
    }
}

OpenHurricane::real OpenHurricane::chemistrySource::tc(const real p, const real T, const realArray &c) {
    reactions().updateG0(T);
    real ciSum = Zero;

    for (integer i = 0; i < nsp(); ++i) {
        ciSum += c[i];
    }
    real ts = Zero;
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        real qf = 0.0;
        real qr = 0.0;
        //qfrj(p, T, c, ri, qf, qr);
        reactions()[ri].q(p, T, c, qf, qr);
        const auto &bc = reactions()[ri].backwardCoeffs();
        for (integer i = 0; i < bc.size(); ++i) {
            integer index = bc[i].index_;
            real stoic = bc[i].stoichCoeff_;
            ts += stoic * qf;
        }
    }
    ts = real(reactions().size()) * ciSum / max(ts, veryTiny);
    if (isnan(ts) || isinf(ts)) {
        return 200;
    }
    return min(ts, real(200));
}

//void OpenHurricane::chemistrySource::DyDt(const real t, const realArray& y, realArray& dydt)
//{
//	dydt = Zero;
//	auto& dcdt = dydt;
//	if (reactionFlowType_ == ReactionFlowType::ConstantPressure)
//	{
//		const real T = y[nsp()];
//		p_ = p0_;
//		for (integer i = 0; i < nsp_; ++i)
//		{
//			yiTmp_[i] = y[i];
//		}
//		auto rho = therm().rho(p_, T, yiTmp_);
//		yiToCi(yiTmp_, rho, ci_);
//		omega(p_, T, ci_, dcdt);
//		real cpmix = Zero;
//		real DTdt = Zero;
//		for (integer i = 0; i < nsp(); i++)
//		{
//			cpmix += yiTmp_[i] * therm()[i].cp_p(p_, T);
//			DTdt += dcdt[i] * therm()[i].Ha_p(p_, T);
//			dydt[i] *= (therm()[i].Wi() / rho);
//		}
//		DTdt /= (rho * cpmix);
//		dcdt[nsp_] = -DTdt;
//		//dcdt[nsp() + 1] = 0.0;
//	}
//	else if (reactionFlowType_ == ReactionFlowType::IsolatedSystem)
//	{
//		errorAbortStr("Not implemented", HUR_FUNCTION);
//	}
//	else if (reactionFlowType_ == ReactionFlowType::ConstantVolume)
//	{
//		const real T = y[nsp_];
//		//reactions().updateG0(T);
//		rho_ = rho0_;
//		for (integer i = 0; i < nsp_; ++i)
//		{
//			yiTmp_[i] = y[i];
//			ci_[i] = rho_ * y[i] / species().W(i);
//		}
//		p_ = therm().p(rho_, T, yiTmp_);
//
//		omega(p_, T, ci_, dcdt);
//
//		real cvmix = Zero;
//		real DTdt = Zero;
//		for (integer i = 0; i < nsp_; i++)
//		{
//			cvmix += yiTmp_[i] * therm()[i].cv_p(p_, T);
//			DTdt += dcdt[i] * therm()[i].ea(rho_, p_, T) * species().W(i);
//			dydt[i] *= (therm()[i].Wi() / rho_);
//		}
//		DTdt /= (cvmix * rho_);
//		dcdt[nsp_] = -DTdt;
//	}
//}

void OpenHurricane::chemistrySource::DyDt(const real t, const realArray &c, realArray &dcdt) {
    dcdt = Zero;
    if (reactionFlowType_ == ReactionFlowType::ConstantPressure) {
        // Temperature
        const real T = c[nsp()];
        p_ = p0_; // The pressure is constant
        for (integer i = 0; i < nsp(); ++i) {
            ci_[i] = max(c[i], real(0));
        }
        omega(p_, T, ci_, dcdt);
        real cCp = Zero;
        real DTdt = Zero;
        for (integer i = 0; i < nsp(); i++) {

            cCp += c[i] * therm()[i].Cp_p(p_, T);

            DTdt += dcdt[i] * therm()[i].Ha_p(p_, T);
        }
        DTdt /= cCp;
        dcdt[nsp()] = -DTdt;

        //dcdt[nsp_ + 1] = 0.0;
    } else if (reactionFlowType_ == ReactionFlowType::ConstantVolume) {
        // Temperature
        const real T = c[nsp()];
        //reactions().updateG0(T);
        rho_ = rho0_; // The density is constant
        for (integer i = 0; i < nsp(); ++i) {
            yiTmp_[i] = c[i] * therm()[i].Wi() / rho_;
            ci_[i] = c[i];
        }
        p_ = therm().p(rho_, T, yiTmp_);

        omega(p_, T, ci_, dcdt);

        real cvmix = Zero;
        real DTdt = Zero;
        for (integer i = 0; i < nsp(); i++) {
            cvmix += yiTmp_[i] * therm()[i].cv_p(p_, T);
            DTdt += dcdt[i] * therm()[i].ea(rho_, p_, T) * species().W(i);
        }
        DTdt /= (cvmix * rho_);
        dcdt[nsp()] = -DTdt;
    }
}

void OpenHurricane::chemistrySource::DyDtSlected(const real t, const realArray &c, realArray &dcdt,
                                             const bool withTemp) {
    dcdt = Zero;
    if (reactionFlowType_ == ReactionFlowType::ConstantPressure) {
        // Temperature
        const real T = c[nsp()];
        p_ = p0_; // The pressure is constant
        for (integer i = 0; i < nsp(); ++i) {
            ci_[i] = max(c[i], real(0));
        }
        omega(p_, T, ci_, dcdt);
        if (withTemp) {
            real cCp = Zero;
            real DTdt = Zero;
            for (integer i = 0; i < nsp(); i++) {
                cCp += c[i] * therm()[i].Cp_p(p_, T);
                DTdt += dcdt[i] * therm()[i].Ha_p(p_, T);
            }
            DTdt /= cCp;
            dcdt[nsp()] = -DTdt;
        }
        //dcdt[nsp_ + 1] = 0.0;
    } else if (reactionFlowType_ == ReactionFlowType::ConstantVolume) {
        // Temperature
        const real T = c[nsp()];
        //reactions().updateG0(T);
        rho_ = rho0_; // The density is constant
        for (integer i = 0; i < nsp(); ++i) {
            yiTmp_[i] = c[i] * therm()[i].Wi() / rho_;
            ci_[i] = c[i];
        }
        p_ = therm().p(rho_, T, yiTmp_);

        omega(p_, T, ci_, dcdt);
        if (withTemp) {
            real cvmix = Zero;
            real DTdt = Zero;
            for (integer i = 0; i < nsp(); i++) {
                cvmix += yiTmp_[i] * therm()[i].cv_p(p_, T);
                DTdt += dcdt[i] * therm()[i].ea(rho_, p_, T) * species().W(i);
            }
            DTdt /= (cvmix * rho_);
            dcdt[nsp()] = -DTdt;
        }
    }
}

void OpenHurricane::chemistrySource::updateOtherValue(const realArray &c) {}

void OpenHurricane::chemistrySource::correctY(realArray &y) {
    real yisum = Zero;
    for (integer i = 0; i < nsp() - 1; ++i) {
        y[i] = max(real(0.0), y[i]);
        yisum += y[i];
    }
    y[nsp() - 1] = 1.0 - max(real(0.0), min(real(1.0), yisum));
}

void OpenHurricane::chemistrySource::jacobian(const real t, const realArray &c, realArray &dfdt,
                                          realSquareMatrix &dfdy) {
    const real T = c[nsp()];
    reactions().updateG0(T);
    for (integer i = 0; i < nsp_; ++i) {
        ci_[i] = max(c[i], real(0));
    }
    dfdy = Zero;
    if (reactionFlowType_ == ReactionFlowType::ConstantVolume) {
        rho_ = rho0_;
        for (integer i = 0; i < nsp_; ++i) {
            yiTmp_[i] = ci_[i] * therm()[i].Wi() / rho_;
        }
        p_ = therm().p(rho_, T, yiTmp_);
    }
    dfdt = Zero;
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        const auto &rr = reactions()[ri];
        real qfr, qf, qr, kf, kr;
        qfr = rr.q(p_, T, ci_, dfdt, qf, qr, kf, kr);
        rr.DqDci(kf, kr, qfr, p_, T, ci_, dfdy);
        rr.DqDT(kf, kr, qfr, p_, T, ci_, dfdy);
    }
    //omega(p_, T, c, dfdt);
    if (reactionFlowType_ == ReactionFlowType::ConstantPressure) {
        real cCp = Zero;
        real DTdt = Zero;
        real dcCpdT = Zero;
        for (integer i = 0; i < nsp(); i++) {
            cCp += c[i] * therm()[i].Cp_p(p_, T);
            dcCpdT += c[i] * therm()[i].DCp_pDT(p_, T);
            DTdt += dfdt[i] * therm()[i].Ha_p(p_, T);
        }
        DTdt /= cCp;
        dfdt[nsp()] = -DTdt;

        // Calculate d(dTdt)/dci
        for (integer i = 0; i < nsp_; ++i) {
            auto &d2Tdtdci = dfdy(nsp_, i);
            for (integer j = 0; j < nsp_; ++j) {
                const real d2cjdtdci = dfdy(j, i);
                d2Tdtdci -= d2cjdtdci * therm()[j].Ha_p(p_, T);
            }
            d2Tdtdci += therm()[i].Cp_p(p_, T) * DTdt;
            d2Tdtdci /= cCp;
        }

        // Calculate d(dTdt)/dT
        real ddTdtdT = Zero;
        for (integer i = 0; i < nsp_; ++i) {
            ddTdtdT -= (therm()[i].Cp_p(p_, T) * dfdt[i] + therm()[i].Ha_p(p_, T) * dfdy(i, nsp_));
        }
        ddTdtdT += dcCpdT * DTdt;
        ddTdtdT /= cCp;
        dfdy(nsp_, nsp_) = ddTdtdT;
    } else if (reactionFlowType_ == ReactionFlowType::ConstantVolume) {
        real cCv = Zero;
        real DTdt = Zero;
        real dcCvdT = Zero;
        for (integer i = 0; i < nsp(); i++) {
            cCv += c[i] * therm()[i].Cv_p(p_, T);
            dcCvdT += c[i] * therm()[i].DCv_pDT(p_, T);
            DTdt += dfdt[i] * therm()[i].ea(rho_, p_, T) * species().W(i);
        }
        DTdt /= cCv;
        dfdt[nsp_] = -DTdt;

        // Calculate d(dTdt)/dc
        for (integer i = 0; i < nsp_; ++i) {
            auto &d2Tdtdci = dfdy(nsp_, i);
            for (integer j = 0; j < nsp_; ++j) {
                const real d2cjdtdci = dfdy(j, i);
                d2Tdtdci -= d2cjdtdci * therm()[j].ea(rho_, p_, T) * species().W(j);
            }
            d2Tdtdci += therm()[i].Cv_p(p_, T) * DTdt;
            d2Tdtdci /= cCv;
        }

        // Calculate d(dTdt)/dT
        real ddTdtdT = Zero;
        for (integer i = 0; i < nsp_; ++i) {
            ddTdtdT -= (therm()[i].Cv_p(p_, T) * dfdt[i] +
                        therm()[i].ea(rho_, p_, T) * species().W(i) * dfdy(i, nsp_));
        }
        ddTdtdT += dcCvdT * DTdt;
        ddTdtdT /= cCv;
        dfdy(nsp_, nsp_) = ddTdtdT;
    }
}

void OpenHurricane::chemistrySource::jacobianSelected(const real t, const realArray &c, realArray &dfdt,
                                                  realSquareMatrix &dfdy, const bool withTemp) {
    const real T = c[nsp()];
    reactions().updateG0(T);
    for (integer i = 0; i < nsp_; ++i) {
        ci_[i] = max(c[i], real(0));
    }
    dfdy = Zero;
    if (reactionFlowType_ == ReactionFlowType::ConstantVolume) {
        rho_ = rho0_;
        for (integer i = 0; i < nsp_; ++i) {
            yiTmp_[i] = ci_[i] * therm()[i].Wi() / rho_;
        }
        p_ = therm().p(rho_, T, yiTmp_);
    }
    dfdt = Zero;
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        const auto &rr = reactions()[ri];
        real qfr, qf, qr, kf, kr;
        qfr = rr.q(p_, T, ci_, dfdt, qf, qr, kf, kr);
        rr.DqDci(kf, kr, qfr, p_, T, ci_, dfdy);
        if (withTemp) {
            rr.DqDT(kf, kr, qfr, p_, T, ci_, dfdy);
        }
    }
    if (!withTemp) {
        return;
    }
    //omega(p_, T, c, dfdt);
    if (reactionFlowType_ == ReactionFlowType::ConstantPressure) {
        real cCp = Zero;
        real DTdt = Zero;
        real dcCpdT = Zero;
        for (integer i = 0; i < nsp(); i++) {
            cCp += c[i] * therm()[i].Cp_p(p_, T);
            dcCpdT += c[i] * therm()[i].DCp_pDT(p_, T);
            DTdt += dfdt[i] * therm()[i].Ha_p(p_, T);
        }
        DTdt /= cCp;
        dfdt[nsp()] = -DTdt;

        // Calculate d(dTdt)/dci
        for (integer i = 0; i < nsp_; ++i) {
            auto &d2Tdtdci = dfdy(nsp_, i);
            for (integer j = 0; j < nsp_; ++j) {
                const real d2cjdtdci = dfdy(j, i);
                d2Tdtdci -= d2cjdtdci * therm()[j].Ha_p(p_, T);
            }
            d2Tdtdci += therm()[i].Cp_p(p_, T) * DTdt;
            d2Tdtdci /= cCp;
        }

        // Calculate d(dTdt)/dT
        real ddTdtdT = Zero;
        for (integer i = 0; i < nsp_; ++i) {
            ddTdtdT -= (therm()[i].Cp_p(p_, T) * dfdt[i] + therm()[i].Ha_p(p_, T) * dfdy(i, nsp_));
        }
        ddTdtdT += dcCpdT * DTdt;
        ddTdtdT /= cCp;
        dfdy(nsp_, nsp_) = ddTdtdT;
    } else if (reactionFlowType_ == ReactionFlowType::ConstantVolume) {
        real cCv = Zero;
        real DTdt = Zero;
        real dcCvdT = Zero;
        for (integer i = 0; i < nsp(); i++) {
            cCv += c[i] * therm()[i].Cv_p(p_, T);
            dcCvdT += c[i] * therm()[i].DCv_pDT(p_, T);
            DTdt += dfdt[i] * therm()[i].ea(rho_, p_, T) * species().W(i);
        }
        DTdt /= cCv;
        dfdt[nsp_] = -DTdt;

        // Calculate d(dTdt)/dc
        for (integer i = 0; i < nsp_; ++i) {
            auto &d2Tdtdci = dfdy(nsp_, i);
            for (integer j = 0; j < nsp_; ++j) {
                const real d2cjdtdci = dfdy(j, i);
                d2Tdtdci -= d2cjdtdci * therm()[j].ea(rho_, p_, T) * species().W(j);
            }
            d2Tdtdci += therm()[i].Cv_p(p_, T) * DTdt;
            d2Tdtdci /= cCv;
        }

        // Calculate d(dTdt)/dT
        real ddTdtdT = Zero;
        for (integer i = 0; i < nsp_; ++i) {
            ddTdtdT -= (therm()[i].Cv_p(p_, T) * dfdt[i] +
                        therm()[i].ea(rho_, p_, T) * species().W(i) * dfdy(i, nsp_));
        }
        ddTdtdT += dcCvdT * DTdt;
        ddTdtdT /= cCv;
        dfdy(nsp_, nsp_) = ddTdtdT;
    }
}

bool OpenHurricane::chemistrySource::checkNewY(real &dt, const realArray &yOld, const realArray &yNew) {
    if ((yNew.last() - yOld.last()) > TemperatureAlpha_ * yOld.last()) {
       /* if (dt < tiny) {
            errorAbortStr(("Timestep size if too small: " + toString(dt) +
                           ". And initial temperature is " + toString(yOld.last()) +
                           ", while new temperature is " + toString(yNew.last())));
        }*/
        dt *= TemperatureBeta_;

        return false;
    }
    return true;
}

void OpenHurricane::chemistrySource::DomegaDciT(const real rho, const real T, const real p,
                                            const realArray &yi, realArray &omegaj,
                                            realSquareMatrix &dwjdciT) const {
    for (integer i = 0; i < nsp_; ++i) {
        ci_[i] = rho * yi[i] / species().W(i);
    }
    omegaj = Zero;
    dwjdciT = Zero;
    for (integer ri = 0; ri < reactions().size(); ++ri) {
        const auto &rr = reactions()[ri];
        real qfr, qf, qr, kf, kr;
        qfr = rr.q(p_, T, ci_, omegaj, qf, qr, kf, kr);
        rr.DqDci(kf, kr, qfr, p_, T, ci_, dwjdciT);
        rr.DqDT(kf, kr, qfr, p_, T, ci_, dwjdciT);
    }
}

OpenHurricane::realArray OpenHurricane::chemistrySource::ciToYi(const realArray &ci, real &rho) const {
    realArray yi(nsp_, Zero);

    for (integer i = 0; i < nsp_; ++i) {
        yi[i] = max(real(0.0), ci[i]);
    }

    rho = Zero;
    for (integer i = 0; i < nsp_; ++i) {
        rho += yi[i] * species().W(i);
    }

    for (integer i = 0; i < nsp_; ++i) {
        yi[i] = yi[i] * species().W(i) / max(rho, tiny);
    }

    return yi;
}

void OpenHurricane::chemistrySource::ciToYi(const realArray &ci, const real rho, realArray &yi) const {
    for (integer i = 0; i < nsp_; ++i) {
        yi[i] = max(real(0.0), ci[i]);
    }

    /*rho = Zero;
    for (integer i = 0; i < nsp_; ++i)
    {
            rho += yi[i] * species().W(i);
    }*/
    real yiSum = 0;
    for (integer i = 0; i < nsp_; ++i) {
        yi[i] = yi[i] * species().W(i) / max(rho, tiny);
        yiSum += yi[i];
    }
    if (yiSum != 1.0) {
        yi /= yiSum;
    }
}

OpenHurricane::real OpenHurricane::chemistrySource::ciToYi(const realArray &ci, realArray &yi) const {
    for (integer i = 0; i < nsp_; ++i) {
        yi[i] = max(real(0.0), ci[i]);
    }

    real rho = Zero;
    for (integer i = 0; i < nsp_; ++i) {
        rho += yi[i] * species().W(i);
    }
    real yiSum = 0;
    for (integer i = 0; i < nsp_; ++i) {
        yi[i] = yi[i] * species().W(i) / max(rho, tiny);
        yiSum += yi[i];
    }
    if (yiSum != 1.0) {
        yi /= yiSum;
    }
    return rho;
}

OpenHurricane::realArray OpenHurricane::chemistrySource::ciToYi(const realArray &ci) const {
    real rho = Zero;
    return ciToYi(ci, rho);
}

OpenHurricane::realArray OpenHurricane::chemistrySource::yiToCi(const realArray &yi, const real rho) const {
    realArray ci(species().size());

    for (integer i = 0; i < ci.size(); ++i) {
        ci[i] = rho * yi[i] / species().W(i);
        ci[i] = max(real(0.0), ci[i]);
    }
    return ci;
}

void OpenHurricane::chemistrySource::yiToCi(const realArray &yi, const real rho, realArray &ci) const {
    if (ci.size() < nsp_) {
        ci.resize(nsp_);
    }
    for (integer i = 0; i < nsp_; ++i) {
        ci[i] = rho * yi[i] / species().W(i);
        ci[i] = max(real(0.0), ci[i]);
    }
}

hur_nodiscard OpenHurricane::real OpenHurricane::chemistrySource::getTemperature(const realArray &ci,
                                                                         const bool check) const {
    if (check) {
        if (ci.size() < this->nsp()) {
            errorAbortStr(("The size of ci is: " + toString(ci.size()) +
                           " is not equal to th esize of species: " + toString(this->nSpc())));
        }
    }
    real rho = Zero;
    for (integer i = 0; i < this->nsp(); ++i) {
        rho += ci[i] * this->species().W(i);
    }
    auto yi = ciToYi(ci, rho);
    real T = this->T_;
    integer iFlag = 0;
    if (this->isConstPressure()) {
        T = this->therm().THa_p(this->hea0_, this->p_, T, iFlag, yi);
        if (check) {
            if (isnan(T)) {
                LFatal("T is nan");
            }
        }
    } else {
        T = this->therm().TEa_rho(this->hea0_, this->rho_, T, iFlag, yi);
        if (check) {
            if (isnan(T)) {
                LFatal("T is nan");
            }
        }
    }
    return T;
}

OpenHurricane::real OpenHurricane::chemistrySource::solve(const real t, const real dT, const real _p,
                                                  const real _T, const real _rho, realArray &_yi) {
    p0_ = _p;
    rho0_ = _rho;
    p_ = p0_;
    rho_ = rho0_;
    T_ = _T;
    real leftTime = t;
    real subDt = dT;
    realArray ci(nsp_);
    for (integer i = 0; i < nsp_; ++i) {
        ci[i] = rho0_ * _yi[i] / therm()[i].Wi();
    }
    while (leftTime > 0.0) {
        real dt = leftTime;
        this->solve(ci, dt, subDt);
        leftTime -= dt;
    }
    ciToYi(ci, _yi);
    if (reactionFlowType_ == ReactionFlowType::ConstantPressure) {
        rho_ = therm().rho(p_, T_, _yi);
    } else {
        p_ = therm().p(rho_, T_, _yi);
    }
    return subDt;
}

OpenHurricane::real OpenHurricane::chemistrySource::solveTest(const real t, const real dT, real &_p,
                                                      real &_T, real &_rho, realArray &_yi,
                                                      fileOsstream &fos) {
    p0_ = _p;
    rho0_ = _rho;
    p_ = p0_;
    rho_ = rho0_;
    T_ = _T;
    real leftTime = t;
    real subDt = dT;
    realArray ci(nsp_);
    for (integer i = 0; i < nsp_; ++i) {
        ci[i] = rho0_ * _yi[i] / therm()[i].Wi();
    }

#ifdef TEST_PROCESS_TIME
    hrClock calcTime;
    odeCountIter_ = 0;
#endif
    while (leftTime > 0.0) {
        real dt = leftTime;
        this->solve(ci, dt, subDt);
        leftTime -= dt;
#ifdef TEST_PROCESS_TIME
        setODECountIter();
#endif
    }
#ifdef TEST_PROCESS_TIME
    calcTime_ += calcTime.clockTimeIncrement();
#endif
    ciToYi(ci, _yi);
    if (reactionFlowType_ == ReactionFlowType::ConstantPressure) {
        rho_ = therm().rho(p_, T_, _yi);
    } else {
        p_ = therm().p(rho_, T_, _yi);
    }
    _rho = rho_;
    _p = p_;
    _T = T_;
    return subDt;
}

void OpenHurricane::chemistrySource::solve(realArray &yi, real &dT, real &subDT) {}

void OpenHurricane::chemistrySource::calculateSourceTerms(const bool withoutLastSpec) {
    if (!isReacting_) {
        return;
    }
    const auto &T = flows_.T();
    const auto &p = flows_.p();
    const auto &rho = flows_.rho();

    for (integer cellI = 0; cellI < mesh_.nCells(); ++cellI) {
        const real Ti = T[cellI];
        const real pi = p[cellI];
        const real rhoi = rho[cellI];

        // Transfer the mass fraction: yi to molar concentration: ci
        for (integer i = 0; i < nsp_; ++i) {
            ci_[i] = rhoi * yi_[i][cellI] / species_.W(i);
        }

        /*if (withoutLastSpec)
        {
                omegaCoupled(pi, Ti, ci_, dcdt_);
                for (integer isp = 0; isp < nsp_ - 1; ++isp)
                {
                        Ri()[isp][cellI] = dcdt_[isp] * species_.W(isp);
                }
        }
        else
        {
                omega(pi, Ti, ci_, dcdt_);
                for (integer isp = 0; isp < nsp_; ++isp)
                {
                        Ri()[isp][cellI] = dcdt_[isp] * species_.W(isp);
                }
        }*/
        omega(pi, Ti, ci_, dcdt_);
        for (integer isp = 0; isp < nsp_; ++isp) {
            Ri()[isp][cellI] = dcdt_[isp] * species_.W(isp);
        }
    }
}

void OpenHurricane::chemistrySource::calculateSourceTermsImp(const bool withoutLastSpec) {
    if (!isReacting_) {
        return;
    }
    const auto &T = flows_.T();
    const auto &p = flows_.p();
    const auto &rho = flows_.rho();
    realArray yyi(nsp_, Zero);
    realArray diagSource(nsp_, Zero);

    for (integer cellI = 0; cellI < mesh_.nCells(); ++cellI) {
        const real Ti = T[cellI];
        const real pi = p[cellI];
        const real rhoi = rho[cellI];

        // Transfer the mass fraction: yi to molar concentration: ci
        for (integer i = 0; i < nsp_; ++i) {
            ci_[i] = rhoi * yi_[i][cellI] / species_.W(i);
            yyi[i] = yi_[i][cellI];
        }

        if (withoutLastSpec) {
            omegaCoupled(rhoi, pi, Ti, ci_, yyi, dcdt_, diagSource);

            for (integer isp = 0; isp < nsp_; ++isp) {
                Ri()[isp][cellI] = dcdt_[isp] * species_.W(isp);
            }
            for (integer isp = 0; isp < nsp_ - 1; ++isp) {
                yi_[isp].diagSource()[cellI] = diagSource[isp];
            }
        } else {
            omega(rhoi, pi, Ti, ci_, yyi, dcdt_, diagSource);
            for (integer isp = 0; isp < nsp_; ++isp) {
                Ri()[isp][cellI] = dcdt_[isp] * species_.W(isp);
                yi_[isp].diagSource()[cellI] = diagSource[isp];
            }
        }
    }
}

void OpenHurricane::chemistrySource::sourceTermsForCellI(const integer cellI, realSquareMatrix &dwdy) {
    if (!isReacting_) {
        return;
    }
    const auto &T = flows_.T();
    const auto &p = flows_.p();
    const auto &rho = flows_.rho();

    const real Ti = T[cellI];
    const real pi = p[cellI];
    const real rhoi = rho[cellI];

    for (integer i = 0; i < nsp_; ++i) {
        ci_[i] = rhoi * yi_[i][cellI] / species_.W(i);
    }
    dwdy.setZero();
    omegaFullJacobian(pi, Ti, ci_, dcdt_, dwdy);
    for (integer isp = 0; isp < nsp_; ++isp) {
        Ri()[isp][cellI] = dcdt_[isp] * species_.W(isp);
    }
}

OpenHurricane::realArray OpenHurricane::chemistrySource::heatReleaseRate() {
    realArray hr(mesh_.nTotalCells(), Zero);
    for (integer isp = 0; isp < nsp_; ++isp) {
        for (integer cellI = 0; cellI < mesh_.nCells(); ++cellI) {
            hr[cellI] -= Ri()[isp][cellI] * therm()[isp].hc();
        }
    }
    return hr;
}

OpenHurricane::real OpenHurricane::chemistrySource::heatReleaseRate(const real p, const real T,
                                                            const realArray &c) {
    const auto dcdt = omega(p, T, c);
    real hrr = 0;
    for (integer isp = 0; isp < nsp_; ++isp) {
        hrr -= dcdt[isp] * species_.W(isp) * therm()[isp].hc();
    }
    return hrr;
}

OpenHurricane::real OpenHurricane::chemistrySource::solve(const realArray &timeStep, const real dtFactor) {
    //#ifdef TEST_PROCESS_TIME
    //	calcTime_ = 0;
    //#endif
    real dtMin = large;
    if (!isReacting_) {
        return dtMin;
    }
    const auto &T = flows_.T();
    const auto &p = flows_.p();
    const auto &rho = flows_.rho();
    realArray ci0(nsp_);
    realArray ci(nsp_);
    auto &dtI = dtInit(timeStep, dtFactor);

    auto &sorCalcTime = chemistrySource::cellSourceCalTime();

    for (integer cellI = 0; cellI < mesh_.nCells(); ++cellI) {
        const real Ti = T[cellI];
        if (Ti < TemperatureThreshold_) {
            dtI[cellI] = timeStep[cellI] * dtFactor;
            for (integer isp = 0; isp < nsp_; ++isp) {
                Ri()[isp][cellI] = 0;
            }
            sorCalcTime[cellI] = 0;
            continue;
        }
        const real pi = p[cellI];
        const real rhoi = rho[cellI];
        hrClock calcSor;

        for (integer i = 0; i < nsp_; ++i) {
            ci[i] = rhoi * yi_[i][cellI] / species_.W(i);
            ci0[i] = ci[i];
        }
        this->p0_ = this->p_ = pi;
        this->rho0_ = this->rho_ = rhoi;
        T_ = Ti;
        real leftTime = timeStep[cellI] * dtFactor;
        real subDt = min(leftTime, min(dtI[cellI], dtMax_));

#ifdef TEST_PROCESS_TIME
        hrClock myClo;
#endif
        while (leftTime > 0.0) {
            real dt = leftTime;
            this->solve(ci, dt, subDt);
            leftTime -= dt;
        }

#ifdef TEST_PROCESS_TIME
        calcTime_ += myClo.elapsedClockTime();
#endif
        dtI[cellI] = subDt;
        dtMin = min(dtMin, subDt);
        for (integer isp = 0; isp < nsp_; ++isp) {
            Ri()[isp][cellI] =
                (ci[isp] - ci0[isp]) / (timeStep[cellI] * dtFactor) * species_.W(isp);
        }
        sorCalcTime[cellI] = calcSor.elapsedClockTime();
    }
    return dtMin;
}

void OpenHurricane::chemistrySource::solveUpdate(const realArray &dt, const real dtFactor) {
    //#ifdef TEST_PROCESS_TIME
    //	calcTime_ = 0;
    //#endif
    if (!isReacting_) {
        return;
    }
    auto &T = flows_.T();
    const auto &p = flows_.p();
    const auto &rho = flows_.rho();

    auto &dtI = dtInit(dt, dtFactor);
    realArray yi0(nsp());

    auto &sorCalcTime = chemistrySource::cellSourceCalTime();

    for (integer cellI = 0; cellI < mesh_.nCells(); ++cellI) {
        const real Ti = T[cellI];
        if (Ti < TemperatureThreshold_) {
            dtI[cellI] = dt[cellI] * dtFactor;

            sorCalcTime[cellI] = 0;
            continue;
        }
        const real pi = p[cellI];
        const real rhoi = rho[cellI];

        hrClock calcSor;

        for (integer i = 0; i < nsp_; ++i) {
            yi0[i] = yi_[i][cellI];
        }
        this->p0_ = this->p_ = pi;
        this->rho0_ = this->rho_ = rhoi;
        T_ = Ti;
        real leftTime = dt[cellI] * dtFactor;
        real subDt = min(leftTime, min(dtI[cellI], dtMax_));

#ifdef TEST_PROCESS_TIME
        hrClock myClo;
#endif
        dtI[cellI] = this->solve(leftTime, subDt, pi, Ti, rhoi, yi0);
#ifdef TEST_PROCESS_TIME
        calcTime_ += myClo.elapsedClockTime();
#endif
        for (integer isp = 0; isp < nsp_; ++isp) {
            yi_[isp][cellI] = yi0[isp];
        }
        T[cellI] = T_;
        sorCalcTime[cellI] = calcSor.elapsedClockTime();
    }
}

//void OpenHurricane::chemistrySource::solve(
//	cellRealArray& rhoi,
//	cellRealArray& pi,
//	cellRealArray& Ti,
//	PtrList<cellRealArray>& yii,
//	const realArray& dt,
//	realArray& subDt,
//	const real dtFactor
//) {
//	if (!isReacting_) {
//		return;
//	}
//	realArray yi0(nsp());
//	for (integer cellI = 0; cellI < mesh_.nCells(); ++cellI)
//	{
//		auto& Tii = Ti[cellI];
//		if (Tii < TemperatureThreshold_)
//		{
//			continue;
//		}
//		auto& pii = pi[cellI];
//		auto& rhoii = rhoi[cellI];
//		for (integer i = 0; i < nsp(); ++i) {
//			yi0[i] = yii[i][cellI];
//		}
//		const auto dtt = dt[cellI] * dtFactor;
//		if (dtt == 0) {
//			errorAbortStr(
//				"The time for reacting is 0 : dt = " + toString(dt[cellI]) + ", dtFactor = " + toString(dtFactor),
//				HUR_FUNCTION
//			);
//		}
//		subDt[cellI] = min(dtt, subDt[cellI]);
//		if (subDt[cellI] == 0) {
//			subDt[cellI] = dtt;
//		}
//
//		subDt[cellI] = this->solve(
//			dtt,
//			subDt[cellI],
//			pii, Tii, rhoii, yi0
//		);
//
//		for (integer i = 0; i < nsp(); ++i) {
//			yii[i][cellI] = yi0[i];
//		}
//	}
//}

void OpenHurricane::chemistrySource::setCellLoadWeights() {}