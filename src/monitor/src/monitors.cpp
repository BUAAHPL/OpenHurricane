/*!
 * \file monitors.cpp
 * \brief Main subroutines for monitors.
 * \author Rao Sihang
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

#include "monitors.hpp"

void OpenHurricane::monitors::removeDuplicate(stringList &nameList) const {
    if (nameList.size() == 0) {
        return;
    }
    integer m = 0;
    stringList tmpList(nameList.size());

    for (integer i = 0; i < nameList.size(); ++i) {
        integer j = 0;
        for (; j < m; ++j) {
            if (nameList[i] == tmpList[j]) {
                break;
            }
        }

        if (j == m) {
            tmpList[m] = nameList[i];
            m++;
        }
    }
    tmpList.resize(m);
    nameList.transfer(tmpList);
}

OpenHurricane::monitors::monitors(const iteration &iter, const runtimeMesh &mesh)
    : iter_(iter), mesh_(mesh), monitorList_() {
    if (iter.cont().subController("iteration").found("monitor")) {
        const auto &monitorCont = iter.cont().subController("iteration").subController("monitor");
        setMonitorList(monitorCont);
    }
}

void OpenHurricane::monitors::setMonitorList(const controller &monitorCont) {
    if (monitorCont.found("monitorList")) {
        std::string rnl = monitorCont.findText("monitorList");
        replaceAllMarks(rnl, "\n", " ");
        stringList monitorList;
        if (!rnl.empty()) {
            size_t pos = 0;
            stdStringList rll;
            split(rnl, rll, ",");
            monitorList.resize(rll.size());
            for (integer i = 0; i < rll.size(); ++i) {
                monitorList[i] = trimCopy(rll[i]);
            }
        }
        removeDuplicate(monitorList);
        if (monitorList_.size() == 0) {
            for (integer i = 0; i < monitorList.size(); ++i) {
                if (monitorCont.found(monitorList[i])) {
                    const auto &icont = monitorCont.subController(monitorList[i]);
                    monitorList_.append(
                        monitor::creator(iter_, mesh_, icont, monitorList[i]).release());
                }
            }
        } else {
            for (integer i = 0; i < monitorList.size(); ++i) {
                if (monitorCont.found(monitorList[i])) {
                    const auto &icont = monitorCont.subController(monitorList[i]);
                    integer m = 0;
                    bool found = false;
                    for (integer j = 0; j < monitorList_.size(); ++j) {
                        if (monitorList_[j].monitorName() == monitorList[i]) {
                            found = true;
                            m = j;
                            break;
                        }
                    }
                    if (found) {
                        monitorList_.set(
                            m, monitor::creator(iter_, mesh_, icont, monitorList[i]).release());
                    } else {
                        monitorList_.append(
                            monitor::creator(iter_, mesh_, icont, monitorList[i]).release());
                    }
                }
            }
        }
    }
}

void OpenHurricane::monitors::monitoring() const {
    if (monitorList_.size() != 0) {
        for (integer i = 0; i < monitorList_.size(); ++i) {
            monitorList_[i].monitoring();
        }
    }
}

void OpenHurricane::monitors::subMonitoring() const {
    if (monitorList_.size() != 0) {
        for (integer i = 0; i < monitorList_.size(); ++i) {
            monitorList_[i].subMonitoring();
        }
    }
}
