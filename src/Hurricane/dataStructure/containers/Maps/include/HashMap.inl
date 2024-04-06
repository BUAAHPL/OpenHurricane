/*!
 * \file HashMap.inl
 * \brief In-Line subroutines of the <i>HashMap.hpp</i> file.
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
#include "HashMap.hpp"
template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq,
                              usePrimeNumber>::HashBucketNode::HashBucketNode(const keyType &key,
                                                                              const Type &ele,
                                                                              HashBucketNode *next)
    : key_(key), elemData_(ele), next_(next) {}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq,
                              usePrimeNumber>::HashBucketNode::HashBucketNode(keyType &&key,
                                                                              const Type &ele,
                                                                              HashBucketNode *next)
    : key_(std::move(key)), elemData_(ele), next_(next) {}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq,
                              usePrimeNumber>::HashBucketNode::HashBucketNode(const keyType &key,
                                                                              Type &&ele,
                                                                              HashBucketNode *next)
    : key_(key), elemData_(std::move(ele)), next_(next) {}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq,
                              usePrimeNumber>::HashBucketNode::HashBucketNode(keyType &&key,
                                                                              Type &&ele,
                                                                              HashBucketNode *next)
    : key_(std::move(key)), elemData_(std::move(ele)), next_(next) {}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::HashMap()
    : HashMapBase(), elementCount_(0), bucketCount_(0), map_(nullptr), maxLoadFactor_(0.8f) {}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::HashMap(
    const HashMap &right)
    : HashMapBase(right), elementCount_(0), bucketCount_(right.bucketCount_), map_(nullptr),
      maxLoadFactor_(right.maxLoadFactor_) {
    if (bucketCount_ != 0) {
        map_ = new HashBucketNode *[bucketCount_];

        for (size_type i = 0; i < bucketCount_; ++i) {
            map_[i] = nullptr;
        }

        for (const_iterator iter = right.cbegin(); iter != right.cend(); ++iter) {
            emplace(iter.key(), iter.element());
        }
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::operator=(
    const HashMap &right) {
    if (this != std::addressof(right)) {
        if (bucketCount_ != 0) {
            clearRelease();
        }
        elementCount_ = 0;
        bucketCount_ = right.bucketCount_;
        if (bucketCount_ != 0) {
            map_ = new HashBucketNode *[bucketCount_];

            for (size_type i = 0; i < bucketCount_; ++i) {
                map_[i] = nullptr;
            }

            for (const_iterator iter = right.cbegin(); iter != right.cend(); ++iter) {
                emplace(iter.key(), iter.element());
            }
        }
        maxLoadFactor_ = right.maxLoadFactor_;
    }
    return *this;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::HashMap(
    HashMap &&right) noexcept
    : HashMapBase(right), elementCount_(right.elementCount_), bucketCount_(right.bucketCount_),
      map_(right.map_), maxLoadFactor_(right.maxLoadFactor_) {
    right.elementCount_ = 0;
    right.bucketCount_ = 0;
    right.map_ = nullptr;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::operator=(
    HashMap &&right) noexcept {
    if (this != std::addressof(right)) {
        if (bucketCount_ != 0) {
            clearRelease();
        }
        elementCount_ = right.elementCount_;
        right.elementCount_ = 0;

        bucketCount_ = right.bucketCount_;
        right.bucketCount_ = 0;

        map_ = right.map_;
        right.map_ = nullptr;

        maxLoadFactor_ = right.maxLoadFactor_;
    }
    return *this;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::HashMap(
    const size_type nbuckets, const float maxLFct)
    : HashMapBase(), elementCount_(0), bucketCount_(0), map_(nullptr), maxLoadFactor_(maxLFct) {
    rehash(nbuckets);
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::HashMap(
    std::initializer_list<std::pair<keyType, Type>> IList)
    : HashMapBase(), elementCount_(0), bucketCount_(0), map_(nullptr), maxLoadFactor_(0.8f) {
    rehash(IList.size());

    for (const std::pair<keyType, Type> &e : IList) {
        emplace(e.first, e.second);
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::operator=(
    std::initializer_list<std::pair<keyType, Type>> IList) {
    if (bucketCount_ != 0) {
        clearRelease();
    }
    elementCount_ = 0;
    bucketCount_ = IList.size();
    if (bucketCount_ != 0) {
        map_ = new HashBucketNode *[bucketCount_];

        for (size_type i = 0; i < bucketCount_; ++i) {
            map_[i] = nullptr;
        }

        for (const std::pair<keyType, Type> &e : IList) {
            emplace(e.first, e.second);
        }
    }

    return *this;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::~HashMap() noexcept {
    if (map_ != nullptr) {
        clear();
        delete[] map_;
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::empty() const noexcept {
    return elementCount_ == 0;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline void OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::clear() noexcept {
    if (elementCount_ == 0) {
        return;
    }
    for (size_type i = 0; i < bucketCount_; ++i) {
        if (map_[i] != nullptr) {
            HashBucketNode *curNode = map_[i];
            HashBucketNode *nextNode = curNode->next_;
            while (nextNode != nullptr) {
                delete curNode;
                curNode = nextNode;
                nextNode = curNode->next_;
            }
            delete curNode;
            map_[i] = nullptr;
        }
    }
    elementCount_ = 0;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline void
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::clearRelease() noexcept {
    clear();
    if (map_ != nullptr) {
        delete[] map_;
        map_ = nullptr;
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::hasher
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::hash_function() const {
    return hasher();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::key_equal
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::key_eq() const {
    return key_equal();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline void OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::rehash(
    const size_type nbuckets) {
    auto newSize = mapSize(nbuckets, isUsingPrime);
    if (nbuckets < 0 || newSize < 0) {
        LFatal("Invalid hash bucket size");
    }
    if (newSize == bucketCount_) {
        return;
    }
    if (newSize == 0) {
        clearRelease();
        return;
    }

    HashBucketNode **newMap = new HashBucketNode *[newSize];
    size_type newEleCount = 0;
    for (size_type i = 0; i < newSize; ++i) {
        newMap[i] = nullptr;
    }
    if (bucketCount_ != 0 && elementCount_ != 0) {
        for (iterator iter = this->begin(), last = this->end(); iter != last;) {
            newEleCount++;
            size_type mapId = mapKeyIndex(iter.key(), newSize);
            HashBucketNode *oldNextNode = iter.curNode_->next_;
            if (oldNextNode != nullptr) {
                if (newMap[mapId] == nullptr) {
                    newMap[mapId] = iter.curNode_;
                    newMap[mapId]->next_ = nullptr;
                } else {
                    HashBucketNode *preNode = newMap[mapId];
                    HashBucketNode *curNode = preNode->next_;
                    while (curNode != nullptr) {
                        preNode = curNode;
                        curNode = curNode->next_;
                    }
                    preNode->next_ = iter.curNode_;
                    preNode->next_->next_ = nullptr;
                }
                iter.curNode_ = oldNextNode;
            } else {
                if (newMap[mapId] == nullptr) {
                    newMap[mapId] = iter.curNode_;
                    newMap[mapId]->next_ = nullptr;
                } else {
                    HashBucketNode *preNode = newMap[mapId];
                    HashBucketNode *curNode = preNode->next_;
                    while (curNode != nullptr) {
                        preNode = curNode;
                        curNode = curNode->next_;
                    }
                    preNode->next_ = iter.curNode_;
                    preNode->next_->next_ = nullptr;
                }
                ++iter;
            }
        }
        for (size_type i = 0; i < bucketCount_; ++i) {
            map_[i] = nullptr;
        }
    }
    HashBucketNode **oldMap = map_;
    map_ = newMap;
    bucketCount_ = newSize;
    elementCount_ = newEleCount;
    delete[] oldMap;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline float
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::load_factor() const noexcept {
    return static_cast<float>(size()) / static_cast<float>(bucket_count());
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline float
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::loadFactor() const noexcept {
    return load_factor();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline float
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::max_load_factor()
    const noexcept {
    return maxLoadFactor_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline float
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::maxLoadFactor()
    const noexcept {
    return max_load_factor();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline void
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::max_load_factor(
    const float factor) const noexcept {
    if (factor <= 0 && std::isnan(factor)) {
        LFatal("Invalid hash load factor");
    }
    maxLoadFactor_ = factor;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline void
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::maxLoadFactor(
    const float factor) const noexcept {
    max_load_factor(factor);
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline void OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::swap(
    HashMap &right) noexcept {
    if (this != std::addressof(right)) {
        auto tmpEleCount = elementCount_;
        elementCount_ = right.elementCount_;
        right.elementCount_ = tmpEleCount;

        auto tmpBucketCount = bucketCount_;
        bucketCount_ = right.bucketCount_;
        right.bucketCount_ = tmpBucketCount;

        auto tmpMap = map_;
        map_ = right.map_;
        right.map_ = tmpMap;

        auto tmpMaxLF = maxLoadFactor_;
        maxLoadFactor_ = right.maxLoadFactor_;
        right.maxLoadFactor_ = tmpMaxLF;
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::iterator(
    HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> *map)
    : mapPtr_(map), mapIndex_(0), curNode_(nullptr) {
    if (mapPtr_ != nullptr) {
        if (mapPtr_->elementCount_ != 0) {
            for (; mapIndex_ < mapPtr_->bucketCount_; ++mapIndex_) {
                if (mapPtr_->map_[mapIndex_] != nullptr) {
                    curNode_ = mapPtr_->map_[mapIndex_];
                    break;
                }
            }
            if (mapIndex_ >= mapPtr_->bucketCount_) {
                mapIndex_ = 0;
            }
        }
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::iterator(
    HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> *map, const size_type mapIndex,
    HashBucketNode *curNode)
    : mapPtr_(map), mapIndex_(mapIndex), curNode_(curNode) {}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::iterator()
    : mapPtr_(nullptr), mapIndex_(0), curNode_(nullptr) {}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline void
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::operator=(
    const iterator &iter) {
    mapPtr_ = iter.mapPtr_;
    mapIndex_ = iter.mapIndex_;
    curNode_ = iter.curNode_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::operator==(
    const iterator &iter) const {
    if (curNode_ == iter.curNode_) {
        return true;
    }
    return false;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::operator!=(
    const iterator &iter) const {
    if (curNode_ == iter.curNode_) {
        return false;
    }
    return true;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::operator++() {
    if (curNode_ != nullptr) {
        if (curNode_->next_ != nullptr) {
            curNode_ = curNode_->next_;
            return *this;
        }
    }
    mapIndex_++;
    for (; mapIndex_ < mapPtr_->bucketCount_; ++mapIndex_) {
        if (mapPtr_->map_[mapIndex_] != nullptr) {
            curNode_ = mapPtr_->map_[mapIndex_];
            break;
        }
    }
    if (mapIndex_ >= mapPtr_->bucketCount_) {
        mapIndex_ = 0;
        curNode_ = nullptr;
    }
    return *this;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::operator++(int) {
    iterator tmp = *this;
    this->operator++();
    return tmp;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::operator*() {
    return curNode_->elemData_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::operator*() const {
    return curNode_->elemData_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq,
                                                     usePrimeNumber>::HashBucketNode &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::node() noexcept {
    return *curNode_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const keyType &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::key()
    const noexcept {
    return curNode_->key_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const keyType &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::first()
    const noexcept {
    return curNode_->key_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::element()
    const noexcept {
    return curNode_->elemData_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::element() noexcept {
    return curNode_->elemData_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::second()
    const noexcept {
    return curNode_->elemData_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator::second() noexcept {
    return curNode_->elemData_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::begin() noexcept {
    return iterator(this);
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::end() noexcept {
    return iterator();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::
    const_iterator(const HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> *map)
    : mapPtr_(map), mapIndex_(0), curNode_(nullptr) {
    if (mapPtr_ != nullptr) {
        if (mapPtr_->elementCount_ != 0) {
            for (; mapIndex_ < mapPtr_->bucketCount_; ++mapIndex_) {
                if (mapPtr_->map_[mapIndex_] != nullptr) {
                    curNode_ = mapPtr_->map_[mapIndex_];
                    break;
                }
            }
            if (mapIndex_ >= mapPtr_->bucketCount_) {
                mapIndex_ = 0;
            }
        }
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::
    const_iterator(const HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> *map,
                   const size_type mapIndex, const HashBucketNode *curNode)
    : mapPtr_(map), mapIndex_(mapIndex), curNode_(curNode) {}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq,
                              usePrimeNumber>::const_iterator::const_iterator()
    : mapPtr_(nullptr), mapIndex_(0), curNode_(nullptr) {}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline OpenHurricane::HashMap<keyType, Type, Hasher, keyEq,
                              usePrimeNumber>::const_iterator::const_iterator(const iterator &iter)
    : mapPtr_(iter.mapPtr_), mapIndex_(iter.mapIndex_), curNode_(iter.curNode_) {}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline void
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::operator=(
    const const_iterator &iter) {
    mapPtr_ = iter.mapPtr_;
    mapIndex_ = iter.mapIndex_;
    curNode_ = iter.curNode_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::operator==(
    const const_iterator &iter) const {
    if (curNode_ == iter.curNode_) {
        return true;
    }
    return false;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::operator!=(
    const const_iterator &iter) const {
    if (curNode_ == iter.curNode_) {
        return false;
    }
    return true;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq,
                                       usePrimeNumber>::const_iterator &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::operator++() {
    if (curNode_ != nullptr) {
        if (curNode_->next_ != nullptr) {
            curNode_ = curNode_->next_;
            return *this;
        }
    }
    mapIndex_++;
    for (; mapIndex_ < mapPtr_->bucketCount_; ++mapIndex_) {
        if (mapPtr_->map_[mapIndex_] != nullptr) {
            curNode_ = mapPtr_->map_[mapIndex_];
            break;
        }
    }
    if (mapIndex_ >= mapPtr_->bucketCount_) {
        mapIndex_ = 0;
        curNode_ = nullptr;
    }
    return *this;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::operator++(
    int) {
    const_iterator tmp = *this;
    this->operator++();
    return tmp;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::operator*()
    const {
    return curNode_->elemData_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq,
                                                           usePrimeNumber>::HashBucketNode &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::node()
    const noexcept {
    return *curNode_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const keyType &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::key()
    const noexcept {
    return curNode_->key_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const keyType &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::first()
    const noexcept {
    return curNode_->key_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::element()
    const noexcept {
    return curNode_->elemData_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator::second()
    const noexcept {
    return curNode_->elemData_;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::begin() const noexcept {
    return const_iterator(this);
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::end() const noexcept {
    return const_iterator();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::cbegin() const noexcept {
    return const_iterator(this);
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::cend() const noexcept {
    return const_iterator();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::find(
        const keyType &keyval) {
    if (elementCount_ != 0) {
        size_type mapIndex = mapKeyIndex(keyval);
        if (map_[mapIndex] != nullptr) {
            HashBucketNode *curNode = map_[mapIndex];
            while (curNode != nullptr) {
                if (keyEq()(curNode->key_, keyval)) {
                    return iterator(this, mapIndex, curNode);
                }
                curNode = curNode->next_;
            }
        }
    }
    return iterator();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::find(
        const keyType &keyval) const {
    if (elementCount_ != 0) {
        size_type mapIndex = mapKeyIndex(keyval);
        if (map_[mapIndex] != nullptr) {
            const HashBucketNode *curNode = map_[mapIndex];
            while (curNode != nullptr) {
                if (keyEq()(curNode->key_, keyval)) {
                    return const_iterator(this, mapIndex, curNode);
                }
                curNode = curNode->next_;
            }
        }
    }
    return const_iterator();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::find(
        const keyType &keyval, const size_type keyIndex) {
    if (elementCount_ != 0) {
        size_type mapIndex = keyIndex;
        if (map_[mapIndex] != nullptr) {
            HashBucketNode *curNode = map_[mapIndex];
            while (curNode != nullptr) {
                if (keyEq()(curNode->key_, keyval)) {
                    return iterator(this, mapIndex, curNode);
                }
                curNode = curNode->next_;
            }
        }
    }
    return iterator();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::const_iterator
    OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::find(
        const keyType &keyval, const size_type keyIndex) const {
    if (elementCount_ != 0) {
        size_type mapIndex = keyIndex;
        if (map_[mapIndex] != nullptr) {
            const HashBucketNode *curNode = map_[mapIndex];
            while (curNode != nullptr) {
                if (keyEq()(curNode->key_, keyval)) {
                    return const_iterator(this, mapIndex, curNode);
                }
                curNode = curNode->next_;
            }
        }
    }
    return const_iterator();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::found(const keyType &keyval) {
    if (elementCount_ != 0) {
        size_type mapIndex = mapKeyIndex(keyval);
        if (map_[mapIndex] != nullptr) {
            HashBucketNode *curNode = map_[mapIndex];
            while (curNode != nullptr) {
                if (keyEq()(curNode->key_, keyval)) {
                    return true;
                }
                curNode = curNode->next_;
            }
        }
    }
    return false;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::found(
    const keyType &keyval) const {
    if (elementCount_ != 0) {
        size_type mapIndex = mapKeyIndex(keyval);
        if (map_[mapIndex] != nullptr) {
            HashBucketNode *curNode = map_[mapIndex];
            while (curNode != nullptr) {
                if (keyEq()(curNode->key_, keyval)) {
                    return true;
                }
                curNode = curNode->next_;
            }
        }
    }
    return false;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::found(const keyType &keyval,
                                                                            size_type &keyIndex) {
    if (elementCount_ != 0) {
        size_type mapIndex = mapKeyIndex(keyval);
        keyIndex = mapIndex;
        if (map_[mapIndex] != nullptr) {
            HashBucketNode *curNode = map_[mapIndex];
            while (curNode != nullptr) {
                if (keyEq()(curNode->key_, keyval)) {
                    return true;
                }
                curNode = curNode->next_;
            }
        }
    }
    return false;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::found(
    const keyType &keyval, size_type &keyIndex) const {
    if (elementCount_ != 0) {
        size_type mapIndex = mapKeyIndex(keyval);
        keyIndex = mapIndex;
        if (map_[mapIndex] != nullptr) {
            HashBucketNode *curNode = map_[mapIndex];
            while (curNode != nullptr) {
                if (keyEq()(curNode->key_, keyval)) {
                    return true;
                }
                curNode = curNode->next_;
            }
        }
    }
    return false;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::contains(
    const keyType &keyval) const {
    return found(keyval);
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline std::pair<
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator, bool>
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::emplace(const keyType &keyval,
                                                                              const Type &val) {
    if (bucketCount_ == 0) {
        rehash(2);
    }

    const size_type mapId = mapKeyIndex(keyval);
    //HashBucketNode* prevNode = nullptr;
    bool isExist = false;
    HashBucketNode *curNode = nullptr;

    for (curNode = map_[mapId]; curNode != nullptr; curNode = curNode->next_) {
        if (keyEq()(curNode->key_, keyval)) {
            isExist = true;
            break;
        }
        //prevNode = curNode;
    }

    if (isExist) {
        return std::pair<iterator, bool>(iterator(this, mapId, curNode), false);
    } else {
        if (map_[mapId] == nullptr) {
            map_[mapId] = new HashBucketNode(keyval, val, nullptr);
        } else {
            map_[mapId] = new HashBucketNode(keyval, val, map_[mapId]);
        }
        elementCount_++;
        autoReHash();

        return std::pair<iterator, bool>(iterator(this, mapId, map_[mapId]), true);
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline std::pair<
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator, bool>
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::emplace(keyType &&keyval,
                                                                              const Type &val) {
    if (bucketCount_ == 0) {
        rehash(2);
    }

    const size_type mapId = mapKeyIndex(keyval);
    //HashBucketNode* prevNode = nullptr;
    bool isExist = false;
    HashBucketNode *curNode = nullptr;

    for (curNode = map_[mapId]; curNode != nullptr; curNode = curNode->next_) {
        if (keyEq()(curNode->key_, keyval)) {
            isExist = true;
            break;
        }
        //prevNode = curNode;
    }

    if (isExist) {
        return std::pair<iterator, bool>(iterator(this, mapId, curNode), false);
    } else {
        if (map_[mapId] == nullptr) {
            map_[mapId] = new HashBucketNode(std::move(keyval), val, nullptr);
        } else {
            map_[mapId] = new HashBucketNode(std::move(keyval), val, map_[mapId]);
        }
        elementCount_++;
        autoReHash();

        return std::pair<iterator, bool>(iterator(this, mapId, map_[mapId]), true);
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline std::pair<
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator, bool>
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::emplace(const keyType &keyval,
                                                                              Type &&val) {
    if (bucketCount_ == 0) {
        rehash(2);
    }

    const size_type mapId = mapKeyIndex(keyval);
    //HashBucketNode *prevNode = nullptr;
    bool isExist = false;
    HashBucketNode *curNode = nullptr;

    for (curNode = map_[mapId]; curNode != nullptr; curNode = curNode->next_) {
        if (keyEq()(curNode->key_, keyval)) {
            isExist = true;
            break;
        }
        //prevNode = curNode;
    }

    if (isExist) {
        return std::pair<iterator, bool>(iterator(this, mapId, curNode), false);
    } else {
        if (map_[mapId] == nullptr) {
            map_[mapId] = new HashBucketNode(keyval, std::move(val), nullptr);
        } else {
            map_[mapId] = new HashBucketNode(keyval, std::move(val), map_[mapId]);
        }
        elementCount_++;
        autoReHash();

        return std::pair<iterator, bool>(iterator(this, mapId, map_[mapId]), true);
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline std::pair<
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator, bool>
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::emplace(keyType &&keyval,
                                                                              Type &&val) {
    if (bucketCount_ == 0) {
        rehash(2);
    }

    const size_type mapId = mapKeyIndex(keyval);
    //HashBucketNode *prevNode = nullptr;
    bool isExist = false;
    HashBucketNode *curNode = nullptr;

    for (curNode = map_[mapId]; curNode != nullptr; curNode = curNode->next_) {
        if (keyEq()(curNode->key_, keyval)) {
            isExist = true;
            break;
        }
        //prevNode = curNode;
    }

    if (isExist) {
        return std::pair<iterator, bool>(iterator(this, mapId, curNode), false);
    } else {
        if (map_[mapId] == nullptr) {
            map_[mapId] = new HashBucketNode(std::move(keyval), std::move(val), nullptr);
        } else {
            map_[mapId] = new HashBucketNode(std::move(keyval), std::move(val), map_[mapId]);
        }
        elementCount_++;
        autoReHash();

        return std::pair<iterator, bool>(iterator(this, mapId, map_[mapId]), true);
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline std::pair<
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator, bool>
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::insert(const keyType &keyval,
                                                                             const Type &val) {
    return emplace(keyval, val);
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline std::pair<
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator, bool>
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::insert(keyType &&keyval,
                                                                             const Type &val) {
    return emplace(std::move(keyval), val);
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline std::pair<
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator, bool>
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::insert(const keyType &keyval,
                                                                             Type &&val) {
    return emplace(keyval, std::move(val));
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline std::pair<
    typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator, bool>
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::insert(keyType &&keyval,
                                                                             Type &&val) {
    return emplace(std::move(keyval), std::move(val));
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline typename OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::iterator
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::erase(const iterator &iter) {
    if (iter.curNode_ != nullptr) {
        HashBucketNode *curNode = map_[iter.mapIndex_];
        HashBucketNode *prevNode = nullptr;
        while (curNode != nullptr) {
            if (curNode == iter.curNode_) {
                break;
            }
            prevNode = curNode;
            curNode = curNode->next_;
        }
        iterator ite = iter;
        ite.curNode_ = nullptr;
        if (prevNode != nullptr) {
            prevNode->next_ = iter.curNode_->next_;
            elementCount_--;
            HurDelete(curNode);
            if (prevNode->next_ != nullptr) {
                ite.mapIndex_ = iter.mapIndex_;
                ite.curNode_ = prevNode->next_;
            } else {
                ite.mapIndex_ = iter.mapIndex_;
                ite.mapIndex_++;
                for (; ite.mapIndex_ < bucketCount_; ++(ite.mapIndex_)) {
                    if (map_[ite.mapIndex_] != nullptr) {
                        ite.curNode_ = map_[ite.mapIndex_];
                        break;
                    }
                }
                if (ite.mapIndex_ >= bucketCount_) {
                    ite.mapIndex_ = 0;
                    ite.curNode_ = nullptr;
                }
            }
        } else {
            curNode = iter.curNode_;
            map_[iter.mapIndex_] = curNode->next_;
            HurDelete(curNode);
            elementCount_--;
            prevNode = map_[iter.mapIndex_];
            if (prevNode != nullptr) {
                ite.mapIndex_ = iter.mapIndex_;
                ite.curNode_ = prevNode;
            } else {
                ite.mapIndex_ = iter.mapIndex_;
                ite.mapIndex_++;
                for (; ite.mapIndex_ < bucketCount_; ++(ite.mapIndex_)) {
                    if (map_[ite.mapIndex_] != nullptr) {
                        ite.curNode_ = map_[ite.mapIndex_];
                        break;
                    }
                }
                if (ite.mapIndex_ >= bucketCount_) {
                    ite.mapIndex_ = 0;
                    ite.curNode_ = nullptr;
                }
            }
        }

        return ite;
    }
    return end();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline bool
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::erase(const keyType &key) {
    iterator iter = find(key);
    if (iter != end()) {
        erase(iter);
        return true;
    }
    return false;
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline void OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::autoReHash() {
    if (load_factor() > max_load_factor() && bucketCount_ < maxSize_) {
        rehash(2 * bucketCount_);
    }
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::at(const keyType &keyval) {
    iterator iter = find(keyval);
    if (iter == end()) {
        LFatal("Unable to find a key in HashMap");
    }
    return iter.second();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
hur_nodiscard inline const Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::at(
    const keyType &keyval) const {
    const_iterator iter = find(keyval);
    if (iter == cend()) {
        LFatal("Unable to find a key in HashMap");
    }
    return iter.second();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline Type &OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::operator[](
    const keyType &keyval) {
    iterator iter = find(keyval);

    if (iter == end()) {
        std::pair<iterator, bool> riter = emplace(keyval, Type());

        if (!riter->second) {
            LFatal("Unable to find and insert a key in HashMap");
        }
        return riter->first.second();
    }

    return iter.second();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline Type &
OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::operator[](keyType &&keyval) {
    iterator iter = find(keyval);

    if (iter == end()) {
        std::pair<iterator, bool> riter = emplace(std::move(keyval), Type());

        if (!riter->second) {
            LFatal("Unable to find and insert a key in HashMap");
        }
        return riter->first.second();
    }

    return iter.second();
}

template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber>
inline const Type &OpenHurricane::HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>::operator[](
    const keyType &keyval) const {
    const_iterator iter = find(keyval);

    if (iter == end()) {
        LFatal("Unable to find a key in HashMap");
    }

    return iter.second();
}
