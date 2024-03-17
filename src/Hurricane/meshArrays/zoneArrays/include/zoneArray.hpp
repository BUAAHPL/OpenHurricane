/*!
 * \file zoneArray.hpp
 * \brief Headers of the zone Array.
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

#pragma once
#include "Array.hpp"
#include "HurMPI.hpp"
#include "NullRefObj.hpp"
#include "iteration.hpp"
#include "object.hpp"

namespace OpenHurricane {

    /*!\brief The template class of Array based on geometry.*/
    template <class Type, class Zone> class zoneArray : public object, public Array<Type> {
    public:
        using Base = Array<Type>;
        using value_type = typename Base::value_type;
        using reference = typename Base::reference;
        using const_reference = typename Base::const_reference;
        using difference_type = typename Base::difference_type;
        using size_type = typename Base::size_type;

        using elementType = typename Base::elementType;
        using subArray = typename Base::subArray;

    private:
        /*!\brief Const reference to zone.*/
        const Zone &zone_;

        mutable uniquePtr<Base> timeSumPtr_;

    protected:
        void setOutputName(const string &n) {
            if (feature<Type>::nElements_ > 1) {
                string nN;
                auto &outnl = object::outputVarNameL();
                outnl.resize(feature<Type>::nElements_);
                for (int i = 0; i < feature<Type>::nElements_; ++i) {
                    nN += "\"";
                    nN += n;
                    nN += "[";
                    nN += std::to_string(i);
                    nN += "]\"";
                    outnl[i] = n + "[" + std::to_string(i) + "]";
                    if (i < feature<Type>::nElements_ - 1) {
                        nN += ",";
                    }
                }
                object::outputVarName() = nN;
            } else {
                object::outputVarName() = n;
            }
            object::hasSetOutputName() = true;
        }

    public:
        /*!\brief Return a null Array.*/
        hur_nodiscard inline static const zoneArray<Type, Zone> &nullObject() {
            return NullRefObj::nullRef<zoneArray<Type, Zone>>();
        }

        inline zoneArray(const Zone &fZ, const char *_c, const registerTable &_tb)
            : object(string(_c + fZ.name()), _tb), Array<Type>(fZ.size(), Zero), zone_(fZ),
              timeSumPtr_(nullptr) {
            setOutputName(string(_c));
        }

        inline zoneArray(const Zone &fZ, const string &_name, const registerTable &_tb)
            : object(string(_name + fZ.name()), _tb), Array<Type>(fZ.size(), Zero), zone_(fZ),
              timeSumPtr_(nullptr) {
            setOutputName(_name);
        }

        inline zoneArray(const Zone &fZ, const char *_c, const registerTable &_tb, const zero)
            : object(string(_c + fZ.name()), _tb), Array<Type>(fZ.size(), Zero), zone_(fZ),
              timeSumPtr_(nullptr) {
            setOutputName(string(_c));
        }

        inline zoneArray(const Zone &fZ, const string &_name, const registerTable &_tb, const zero)
            : object(string(_name + fZ.name()), _tb), Array<Type>(fZ.size(), Zero), zone_(fZ),
              timeSumPtr_(nullptr) {
            setOutputName(_name);
        }

        inline zoneArray(const Zone &fZ, const char *_c, const registerTable &_tb, const Type &t)
            : object(string(_c + fZ.name()), _tb), Array<Type>(fZ.size(), t), zone_(fZ),
              timeSumPtr_(nullptr) {
            setOutputName(string(_c));
        }

        inline zoneArray(const Zone &fZ, const string &_name, const registerTable &_tb,
                         const Type &t)
            : object(string(_name + fZ.name()), _tb), Array<Type>(fZ.size(), t), zone_(fZ),
              timeSumPtr_(nullptr) {
            setOutputName(_name);
        }

        inline zoneArray(const Zone &fZ, const char *_c, const registerTable &_tb,
                         const Array<Type> &f)
            : object(string(_c + fZ.name()), _tb), Array<Type>(f), zone_(fZ), timeSumPtr_(nullptr) {
            setOutputName(string(_c));
        }

        inline zoneArray(const Zone &fZ, const string &_name, const registerTable &_tb,
                         const Array<Type> &f)
            : object(string(_name + fZ.name()), _tb), Array<Type>(f), zone_(fZ),
              timeSumPtr_(nullptr) {
            setOutputName(_name);
        }

        inline zoneArray(const Zone &fZ, const char *_c, const registerTable &_tb, Array<Type> &&f)
            : object(string(_c + fZ.name()), _tb), Array<Type>(std::move(f)), zone_(fZ),
              timeSumPtr_(nullptr) {
            setOutputName(string(_c));
        }

        inline zoneArray(const Zone &fZ, const string &_name, const registerTable &_tb,
                         Array<Type> &&f)
            : object(string(_name + fZ.name()), _tb), Array<Type>(std::move(f)), zone_(fZ),
              timeSumPtr_(nullptr) {
            setOutputName(_name);
        }

        inline zoneArray(const zoneArray<Type, Zone> &ff)
            : object(ff), Array<Type>(ff), zone_(ff.zone_), timeSumPtr_(nullptr) {}

        inline zoneArray(zoneArray<Type, Zone> &&ff) noexcept
            : object(std::move(ff)), Array<Type>(std::move(ff)), zone_(ff.zone_),
              timeSumPtr_(nullptr) {}

        hur_nodiscard inline uniquePtr<zoneArray<Type, Zone>> clone() const {
            return uniquePtr<zoneArray<Type, Zone>>(new zoneArray<Type, Zone>(*this));
        }

        inline ~zoneArray() noexcept { timeSumPtr_.clear(); }

        /*!\brief Return the const reference to the mesh.*/
        hur_nodiscard inline const Zone &zone() const noexcept { return zone_; }

        /*!\brief Return const access to the Array.*/
        hur_nodiscard inline const Base &array_ref() const noexcept { return *this; }

        /*!\brief Return non-const access to the Array.*/
        hur_nodiscard inline Base &array_ref() noexcept { return *this; }

        hur_nodiscard inline virtual int nElements() const noexcept {
            return feature<Type>::nElements_;
        }

        hur_nodiscard inline const Array<Type> &timeSum() const {
            if (!timeSumPtr_) {
                timeSumPtr_.reset(new Array<Type>(zone_.size(), Zero));
            }
            return *timeSumPtr_;
        }

        hur_nodiscard inline Base &timeSum() {
            if (!timeSumPtr_) {
                timeSumPtr_.reset(new Base(zone_.size(), Zero));
            }
            return *timeSumPtr_;
        }

        inline virtual void calcTimeSumPtr(const real &dt) const {
            if (!timeSumPtr_) {
                timeSumPtr_.reset(new Base(zone_.size(), Zero));
            }
            (*timeSumPtr_) += dt * (*this);
        }

        using Base::average;
        using Base::weightedAverage;

        void clear() noexcept {
            Base::clear();
            object::clear();
            timeSumPtr_.clear();
        }

        // Member operators

        inline zoneArray<Type, Zone> &operator=(const zoneArray<Type, Zone> &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        inline zoneArray<Type, Zone> &operator=(const Base &other) {
            if (this != std::addressof(other)) {
                Base::operator=(other);
            }
            return *this;
        }

        inline zoneArray<Type, Zone> &operator=(zoneArray<Type, Zone> &&other) noexcept {
            object::operator=(std::move(other));
            Base::operator=(std::move(other));
            timeSumPtr_ = std::move(other.timeSumPtr_);
            return *this;
        }
        inline zoneArray<Type, Zone> &operator=(Base &&other) noexcept {
            Base::operator=(std::move(other));
            return *this;
        }

        inline zoneArray<Type, Zone> &operator=(const Type &t) {
            Base::operator=(t);
            return *this;
        }

        inline zoneArray<Type, Zone> &operator=(const zero) {
            Base::operator=(Zero);
            return *this;
        }

        inline void operator+=(const zoneArray<Type, Zone> &other) { Base::operator+=(other); }
        inline void operator+=(const Base &gf) { Base::operator+=(gf); }
        inline void operator+=(const Type &t) { Base::operator+=(t); }

        inline void operator-=(const zoneArray<Type, Zone> &gf) { Base::operator-=(gf); }
        inline void operator-=(const Base &gf) { Base::operator-=(gf); }
        inline void operator-=(const Type &t) { Base::operator-=(t); }

        inline void operator*=(const zoneArray<Type, Zone> &gf) { Base::operator*=(gf); }
        inline void operator*=(const Base &gf) { Base::operator*=(gf); }
        inline void operator*=(const Type &t) { Base::operator*=(t); }

        inline void operator/=(const zoneArray<Type, Zone> &gf) { Base::operator/=(gf); }
        inline void operator/=(const Base &gf) { Base::operator/=(gf); }
        inline void operator/=(const Type &t) { Base::operator/=(t); }

        // Write

        /*!
         * \brief Write object to output file.
         * Note: it is a virtual function, and should be rewritten in derived class.
         */
        inline virtual void writeOutput(fileOsstream &fos) const {
            // Only write the internale filed's value.
            Base::writeToStream(fos, this->size());
        }
        inline virtual void writeOutput(fileOsstream &fos, const integer fzid) const {}

        /*!
         * \brief Write object minimum and maximum value to output file.
         * Note: it is a virtual function, and should be rewritten in derived class.
         */
        inline virtual void writeMinMaxOutput(fileOsstream &fos) const {}
        inline virtual void writeMinMaxOutput(fileOsstream &fos, const integer fzid) const {}
    };

    // Global functions
    template <class Type, class Zone>
    inline void component(zoneArray<typename zoneArray<Type, Zone>::elementType, Zone> &comp,
                          const zoneArray<Type, Zone> gf, const int d) {
        component(dynamic_cast<typename Array<Type>::elementType &>(comp),
                  dynamic_cast<const Array<Type> &>(gf), d);
    }

    template <class Type, class Zone>
    void T(zoneArray<Type, Zone> &tran, const zoneArray<Type, Zone> gf) {
        T(dynamic_cast<typename Array<Type>::elementType &>(tran),
          dynamic_cast<const Array<Type> &>(gf));
    }

    template <class Type, class Zone>
    hur_nodiscard zoneArray<Type, Zone> operator+(const zoneArray<Type, Zone> &f1,
                                                  const zoneArray<Type, Zone> &f2) {
        zoneArray<Type, Zone> gf(f1.faceZoneS(), f1.name() + "+" + f2.name(), f1.tb());

        for (integer i = 0; i < f1.size(); ++i) {
            gf[i] = f1[i] + f2[i];
        }

        return gf;
    }

    template <class Type, class Zone>
    hur_nodiscard zoneArray<Type, Zone> operator+(zoneArray<Type, Zone> &&f1,
                                                  const zoneArray<Type, Zone> &f2) {
        zoneArray<Type, Zone> tf(std::move(f1));
        for (integer i = 0; i < tf.size(); ++i) {
            tf[i] += f2[i];
        }

        return tf;
    }

    template <class Type, class Zone>
    hur_nodiscard zoneArray<Type, Zone> operator+(const zoneArray<Type, Zone> &f1,
                                                  zoneArray<Type, Zone> &&f2) {
        zoneArray<Type, Zone> tf(std::move(f2));
        for (integer i = 0; i < f1.size(); ++i) {
            tf[i] = f1[i] + tf[i];
        }

        return tf;
    }

    template <class Type, class Zone>
    hur_nodiscard zoneArray<Type, Zone> operator+(zoneArray<Type, Zone> &&f1,
                                                  zoneArray<Type, Zone> &&f2) {
        zoneArray<Type, Zone> tf(std::move(f1));
        for (integer i = 0; i < tf.size(); ++i) {
            tf[i] += f2[i];
        }

        return tf;
    }

    template <class Type, class Zone>
    hur_nodiscard zoneArray<Type, Zone> operator-(const zoneArray<Type, Zone> &f1,
                                                  const zoneArray<Type, Zone> &f2) {
        zoneArray<Type, Zone> gf(f1.faceZoneS(), f1.name() + "-" + f2.name(), f1.tb());

        for (integer i = 0; i < f1.size(); ++i) {
            gf[i] = f1[i] - f2[i];
        }

        return gf;
    }

    template <class Type, class Zone>
    hur_nodiscard zoneArray<Type, Zone> operator-(const zoneArray<Type, Zone> &f1,
                                                  zoneArray<Type, Zone> &&f2) {
        zoneArray<Type, Zone> tf(std::move(f2));
        for (integer i = 0; i < f1.size(); ++i) {
            tf[i] = f1[i] - tf[i];
        }

        return tf;
    }

    template <class Type, class Zone>
    hur_nodiscard zoneArray<Type, Zone> operator-(zoneArray<Type, Zone> &&f1,
                                                  const zoneArray<Type, Zone> &f2) {
        zoneArray<Type, Zone> tf(std::move(f1));
        for (integer i = 0; i < tf.size(); ++i) {
            tf[i] -= f2[i];
        }

        return tf;
    }

    template <class Type, class Zone>
    hur_nodiscard zoneArray<Type, Zone> operator-(zoneArray<Type, Zone> &&f1,
                                                  zoneArray<Type, Zone> &&f2) {
        zoneArray<Type, Zone> tf(std::move(f1));
        for (integer i = 0; i < tf.size(); ++i) {
            tf[i] -= f2[i];
        }

        return tf;
    }

} // namespace OpenHurricane
