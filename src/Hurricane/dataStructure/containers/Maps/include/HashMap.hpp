/*!
 * \file HashMap.hpp
 * \brief Header of HashMap.
 *       The subroutines and functions are in the <i>HashMap.inl</i> file.
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
#include "errorAbort.hpp"
#include "integer.hpp"
#include <cmath>

namespace OpenHurricane {
    class HashMapBase {
    public:
        static constexpr uinteger maxSize_ = uintegerMax / 2;

        hur_nodiscard static uinteger mapSize(const uinteger n, const bool usePrime);
    };
    /**
     * \brief The class of HashMap.
     */
    template <class keyType, class Type, class Hasher, class keyEq, bool usePrimeNumber = false>
    class HashMap : public HashMapBase {
    public:
        using hasher = Hasher;
        using key_type = keyType;
        using mapped_type = Type;
        using key_equal = keyEq;

        using value_type = mapped_type;
        using reference = value_type &;
        using const_reference = const value_type &;
        using size_type = uinteger;

        //using HashMapType = typename HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber>;

        static const bool isUsingPrime = usePrimeNumber;

    private:
        /**
         * \breif Hash bucket node by using single linked list.
         */
        class HashBucketNode {
        public:
            keyType key_;

            Type elemData_;

            HashBucketNode *next_;

            /**
             * \brief Construct from key, element data and the next pointer.
             */
            HashBucketNode(const keyType &key, const Type &ele, HashBucketNode *next);

            /**
             * \brief Construct from key, element data and the next pointer.
             */
            HashBucketNode(keyType &&key, const Type &ele, HashBucketNode *next);

            /**
             * \brief Construct from key, element data and the next pointer.
             */
            HashBucketNode(const keyType &key, Type &&ele, HashBucketNode *next);

            /**
             * \brief Construct from key, element data and the next pointer.
             */
            HashBucketNode(keyType &&key, Type &&ele, HashBucketNode *next);
        };

        /**
         * \brief The number of element data.
         */
        size_type elementCount_;

        /**
         * brief The number of buckets allocated.
         */
        size_type bucketCount_;

        HashBucketNode **map_;

        float maxLoadFactor_;

        hur_nodiscard inline size_type mapKeyIndex(const keyType &key) const {
            if (bucketCount_ == 0) {
                return 0;
            }
            if (usePrimeNumber) {
                return size_type(hasher()(key)) % (bucketCount_);
            } else {
                return size_type(hasher()(key)) & (bucketCount_ - 1);
            }
        }

        hur_nodiscard inline size_type mapKeyIndex(const keyType &key,
                                                   const size_type bucketSize) const {
            if (bucketSize == 0) {
                return 0;
            }
            if (usePrimeNumber) {
                return size_type(hasher()(key)) % (bucketSize);
            } else {
                return size_type(hasher()(key)) & (bucketSize - 1);
            }
        }

        hur_nodiscard inline size_type mapKeyIndexWithHashVal(const size_type ikey) const {
            if (bucketCount_ == 0) {
                return 0;
            }
            if (usePrimeNumber) {
                return ikey % (bucketCount_);
            } else {
                return ikey & (bucketCount_ - 1);
            }
        }

        hur_nodiscard inline size_type countBucketSize(HashBucketNode *curNode) const noexcept {
            if (curNode != nullptr) {
                HashBucketNode *curNode0 = curNode;
                HashBucketNode *nextNode = curNode0->next_;

                size_type count = 1;
                while (nextNode != nullptr) {
                    count++;
                    curNode0 = nextNode;
                    nextNode = curNode0->next_;
                }
                return count;
            }
            return 0;
        }

    public:
        inline HashMap();

        inline HashMap(const HashMap &right);

        inline HashMap &operator=(const HashMap &right);

        inline HashMap(HashMap &&right) noexcept;

        inline HashMap &operator=(HashMap &&right) noexcept;

        inline HashMap(const size_type nbuckets, const float maxLFct = 0.8f);

        inline HashMap(std::initializer_list<std::pair<keyType, Type>> IList);

        inline HashMap &operator=(std::initializer_list<std::pair<keyType, Type>> IList);

        virtual ~HashMap() noexcept;

        /**
         * \brief Return the size of the Hash map, i.e., the number of element data.
         */
        hur_nodiscard inline size_type size() const noexcept { return elementCount_; }

        hur_nodiscard inline size_type max_size() const noexcept { return maxSize_; }

        /**
         * \brief Gets the maximum number of buckets.
         */
        hur_nodiscard inline size_type max_bucket_count() const noexcept { return bucket_count(); }

        /**
         * \brief Return true if the linked list is empty.
         */
        hur_nodiscard inline bool empty() const noexcept;

        /**
         * \brief The number of buckets allocated.
         */
        hur_nodiscard inline size_type capacity() const noexcept { return bucketCount_; }

        /**
         * \brief The number of buckets allocated.
         */
        hur_nodiscard inline size_type bucket_count() const noexcept { return bucketCount_; }

        /**
         * \brief Gets the size of a bucket.
         */
        hur_nodiscard inline size_type bucket_size(const size_type ibucket) const noexcept {
            if (bucketCount_ > 0) {
                return countBucketSize(map_[ibucket]);
            }
            return 0;
        }

        /**
         * \brief Finds the number of elements matching a specified key.
         */
        hur_nodiscard inline size_type count(const keyType &key) const {
            if (bucketCount_ > 0) {
                HashBucketNode *curNode = map_[mapKeyIndex(key)];
                while (curNode != nullptr) {
                    if (keyEq()(curNode->key_, key)) {
                        return 1;
                    }
                    curNode = curNode->next_;
                }
            }
            return 0;
        }

        /**
         * \brief Gets the bucket number for a key value.
         */
        hur_nodiscard inline size_type bucket(const keyType &key) const { return mapKeyIndex(key); }

        void clear() noexcept;

        void clearRelease() noexcept;

        hur_nodiscard inline hasher hash_function() const;

        hur_nodiscard inline key_equal key_eq() const;

        /**
         * \brief Rebuilds the hash table.
         * \param[in] nbuckets - The requested number of buckets.
         */
        void rehash(const size_type nbuckets);

        /**
         * \brief Counts the average elements per bucket.
         */
        hur_nodiscard inline float load_factor() const noexcept;

        /**
         * \brief Counts the average elements per bucket.
         */
        hur_nodiscard inline float loadFactor() const noexcept;

        /**
         * \brief Gets the maximum elements per bucket.
         */
        hur_nodiscard inline float max_load_factor() const noexcept;

        /**
         * \brief Gets the maximum elements per bucket.
         */
        hur_nodiscard inline float maxLoadFactor() const noexcept;

        /**
         * \brief Sets the maximum elements per bucket.
         */
        hur_nodiscard inline void max_load_factor(const float factor) const noexcept;

        /**
         * \brief Sets the maximum elements per bucket.
         */
        hur_nodiscard inline void maxLoadFactor(const float factor) const noexcept;

        void swap(HashMap &right) noexcept;

    public:
        /**
         * \brief The class of iterator in HashMap.
         */
        class iterator {
        public:
            using value_type = Type;
            using size_type = uinteger;

        private:
            friend class HashMap;

            HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> *mapPtr_;

            size_type mapIndex_;

            HashBucketNode *curNode_;

        public:
            inline explicit iterator(HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> *map);
            inline iterator(HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> *map,
                            const size_type size_type, HashBucketNode *curNode);

            inline iterator();

            virtual ~iterator() noexcept {}

            inline void operator=(const iterator &iter);
            hur_nodiscard inline bool operator==(const iterator &) const;
            hur_nodiscard inline bool operator!=(const iterator &) const;

            inline iterator &operator++();
            inline iterator operator++(int);

            hur_nodiscard inline Type &operator*();
            hur_nodiscard inline const Type &operator*() const;

            hur_nodiscard inline HashBucketNode &node() noexcept;

            hur_nodiscard inline const keyType &key() const noexcept;
            hur_nodiscard inline const keyType &first() const noexcept;

            hur_nodiscard inline const Type &element() const noexcept;
            hur_nodiscard inline Type &element() noexcept;

            hur_nodiscard inline const Type &second() const noexcept;
            hur_nodiscard inline Type &second() noexcept;
        };

        /**
         * \brief Returns a iterator that addresses the first element in the range.
         */
        hur_nodiscard inline iterator begin() noexcept;

        /**
         * \brief Returns a iterator that addresses the location just beyond the last element in a range.
         */
        hur_nodiscard inline iterator end() noexcept;

        /**
         * \brief The class of const_iterator in HashMap.
         */
        class const_iterator {
        public:
            using value_type = Type;
            using size_type = uinteger;

        private:
            friend class HashMap;
            friend class iterator;

            const HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> *mapPtr_;

            size_type mapIndex_;

            const HashBucketNode *curNode_;

        public:
            inline explicit const_iterator(
                const HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> *map);
            inline const_iterator(const HashMap<keyType, Type, Hasher, keyEq, usePrimeNumber> *map,
                                  const size_type mapIndex, const HashBucketNode *curNode);

            inline const_iterator();

            inline const_iterator(const iterator &);

            virtual ~const_iterator() noexcept {}

            inline void operator=(const const_iterator &iter);
            hur_nodiscard inline bool operator==(const const_iterator &) const;
            hur_nodiscard inline bool operator!=(const const_iterator &) const;

            inline const_iterator &operator++();
            inline const_iterator operator++(int);

            hur_nodiscard inline const Type &operator*() const;

            hur_nodiscard inline const HashBucketNode &node() const noexcept;
            hur_nodiscard inline const keyType &key() const noexcept;
            hur_nodiscard inline const keyType &first() const noexcept;

            hur_nodiscard inline const Type &element() const noexcept;
            hur_nodiscard inline const Type &second() const noexcept;
        };

        /**
         * \brief Returns a const iterator that addresses the first element in the range.
         */
        hur_nodiscard inline const_iterator begin() const noexcept;

        /**
         * \brief Returns a const iterator that addresses the location just beyond the last element in a range.
         */
        hur_nodiscard inline const_iterator end() const noexcept;

        /**
         * \brief Returns a const iterator that addresses the first element in the range.
         */
        hur_nodiscard inline const_iterator cbegin() const noexcept;

        /**
         * \brief Returns a const iterator that addresses the location just beyond the last element in a range.
         */
        hur_nodiscard inline const_iterator cend() const noexcept;

        friend class iterator;
        friend class const_iterator;

        using local_iterator = iterator;
        using const_local_iterator = const_iterator;

        hur_nodiscard inline local_iterator begin(const size_type ibucket) noexcept {
            return local_iterator(this, ibucket, map_[ibucket]);
        }

        hur_nodiscard inline const_local_iterator begin(const size_type ibucket) const noexcept {
            return const_local_iterator(this, ibucket, map_[ibucket]);
        }

        hur_nodiscard inline local_iterator end(const size_type ibucket) noexcept {
            return local_iterator(this, ibucket, nullptr);
        }

        hur_nodiscard inline const_local_iterator end(const size_type ibucket) const noexcept {
            return const_local_iterator(this, ibucket, nullptr);
        }

        hur_nodiscard inline const_local_iterator cbegin(const size_type ibucket) const noexcept {
            return const_local_iterator(this, ibucket, map_[ibucket]);
        }

        hur_nodiscard inline const_local_iterator cend(const size_type ibucket) const noexcept {
            return const_local_iterator(this, ibucket, nullptr);
        }

    public:
        /**
         * \brief Finds an element that matches a specified key.
         * \param[in] keyval - Key value to search for.
         */
        hur_nodiscard iterator find(const keyType &keyval);

        /**
         * \brief Finds an element that matches a specified key.
         * \param[in] keyval - Key value to search for.
         */
        hur_nodiscard const_iterator find(const keyType &keyval) const;

        /**
         * \brief Finds an element that matches a specified key.
         * \param[in] keyval - Key value to search for.
         * \param[in] keyIndex - The index of this key in the bucket map, i.e., the index of the bucket.
         */
        hur_nodiscard iterator find(const keyType &keyval, const size_type keyIndex);

        /**
         * \brief Finds an element that matches a specified key.
         * \param[in] keyval - Key value to search for.
         * \param[in] keyIndex - The index of this key in the bucket map, i.e., the index of the bucket.
         */
        hur_nodiscard const_iterator find(const keyType &keyval, const size_type keyIndex) const;

        /**
         * \brief To find whether there is an element that matches a specified key.
         * \param[in] keyval - Key value to search for.
         */
        hur_nodiscard bool found(const keyType &keyval);

        /**
         * \brief To find whether there is an element that matches a specified key.
         * \param[in] keyval - Key value to search for.
         */
        hur_nodiscard bool found(const keyType &keyval) const;

        /**
         * \brief To find whether there is an element that matches a specified key.
         * \param[in] keyval - Key value to search for.
         * \param[in] keyIndex - The index of this key in the bucket map, i.e., the index of the bucket.
         */
        hur_nodiscard bool found(const keyType &keyval, size_type &keyIndex);

        /**
         * \brief To find whether there is an element that matches a specified key.
         * \param[in] keyval - Key value to search for.
         * \param[in] keyIndex - The index of this key in the bucket map, i.e., the index of the bucket.
         */
        hur_nodiscard bool found(const keyType &keyval, size_type &keyIndex) const;

        /**
         * \brief Checks if there's an element in the HashMap with the specified key.
         * \param[in] keyval - The key value of the element to look for.
         */
        hur_nodiscard inline bool contains(const keyType &keyval) const;

        /**
         * \brief Inserts an element constructed in place into a HashMap.
         * \param[in] keyval - The key of an element to be inserted into the HashMap.
         * \param[in] val - The value of an element to be inserted into the HashMap.
         * \return Return a pair whose bool component is true if an insertion was made, and false if the HashMap already contained an element with the given key.
         * The iterator component of the return-value pair points to the newly inseted element if the bool component is true, or to the existing element if the bool component is false.
         */
        inline std::pair<iterator, bool> emplace(const keyType &keyval, const Type &val);

        /**
         * \brief Inserts an element constructed in place into a HashMap.
         * \param[in] keyval - The key of an element to be inserted into the HashMap.
         * \param[in] val - The value of an element to be inserted into the HashMap.
         * \return Return a pair whose bool component is true if an insertion was made, and false if the HashMap already contained an element with the given key.
         * The iterator component of the return-value pair points to the newly inseted element if the bool component is true, or to the existing element if the bool component is false.
         */
        inline std::pair<iterator, bool> emplace(keyType &&keyval, const Type &val);

        /**
         * \brief Inserts an element constructed in place into a HashMap.
         * \param[in] keyval - The key of an element to be inserted into the HashMap.
         * \param[in] val - The value of an element to be inserted into the HashMap.
         * \return Return a pair whose bool component is true if an insertion was made, and false if the HashMap already contained an element with the given key.
         * The iterator component of the return-value pair points to the newly inseted element if the bool component is true, or to the existing element if the bool component is false.
         */
        inline std::pair<iterator, bool> emplace(const keyType &keyval, Type &&val);

        /**
         * \brief Inserts an element constructed in place into a HashMap.
         * \param[in] keyval - The key of an element to be inserted into the HashMap.
         * \param[in] val - The value of an element to be inserted into the HashMap.
         * \return Return a pair whose bool component is true if an insertion was made, and false if the HashMap already contained an element with the given key.
         * The iterator component of the return-value pair points to the newly inseted element if the bool component is true, or to the existing element if the bool component is false.
         */
        inline std::pair<iterator, bool> emplace(keyType &&keyval, Type &&val);

        /**
         * \brief Inserts an element constructed in place into a HashMap.
         * \param[in] keyval - The key of an element to be inserted into the HashMap.
         * \param[in] val - The value of an element to be inserted into the HashMap.
         * \return Return a pair whose bool component is true if an insertion was made, and false if the HashMap already contained an element with the given key.
         * The iterator component of the return-value pair points to the newly inseted element if the bool component is true, or to the existing element if the bool component is false.
         */
        inline std::pair<iterator, bool> insert(const keyType &keyval, const Type &val);

        /**
         * \brief Inserts an element constructed in place into a HashMap.
         * \param[in] keyval - The key of an element to be inserted into the HashMap.
         * \param[in] val - The value of an element to be inserted into the HashMap.
         * \return Return a pair whose bool component is true if an insertion was made, and false if the HashMap already contained an element with the given key.
         * The iterator component of the return-value pair points to the newly inseted element if the bool component is true, or to the existing element if the bool component is false.
         */
        inline std::pair<iterator, bool> insert(keyType &&keyval, const Type &val);

        /**
         * \brief Inserts an element constructed in place into a HashMap.
         * \param[in] keyval - The key of an element to be inserted into the HashMap.
         * \param[in] val - The value of an element to be inserted into the HashMap.
         * \return Return a pair whose bool component is true if an insertion was made, and false if the HashMap already contained an element with the given key.
         * The iterator component of the return-value pair points to the newly inseted element if the bool component is true, or to the existing element if the bool component is false.
         */
        inline std::pair<iterator, bool> insert(const keyType &keyval, Type &&val);

        /**
         * \brief Inserts an element constructed in place into a HashMap.
         * \param[in] keyval - The key of an element to be inserted into the HashMap.
         * \param[in] val - The value of an element to be inserted into the HashMap.
         * \return Return a pair whose bool component is true if an insertion was made, and false if the HashMap already contained an element with the given key.
         * The iterator component of the return-value pair points to the newly inseted element if the bool component is true, or to the existing element if the bool component is false.
         */
        inline std::pair<iterator, bool> insert(keyType &&keyval, Type &&val);

        /**
         * \brief Removes an element in a HashMap from the specified position.
         * \param[in] iter - Position of the element to be removed.
         * \return A bidirectional iterator that designates the first element remaining beyond any elements removed,
         *  or an element that is the end of the map if no such element exists.
         */
        iterator erase(const iterator &iter);

        /**
         * \brief Removes element that match a specified key.
         * \param[in] key - The key value of the elements to be removed.
         * \return Return true if the elements have been removed from the HashMap.
         */
        bool erase(const keyType &key);

    protected:
        /**
         * \brief Auto rebuild the map when the load factor is lager than the max load factor.
         */
        inline void autoReHash();

    public:
        /**
         * \brief inds an element in a unordered_map with a specified key value.
         * \param[in] keyval - The key value to find.
         * \return A reference to the data value of the element found.
         * \throw If the argument key value isn't found, then the function throws an error
         */
        hur_nodiscard inline Type &at(const keyType &keyval);

        /**
         * \brief inds an element in a unordered_map with a specified key value.
         * \param[in] keyval - The key value to find.
         * \return A reference to the data value of the element found.
         * \throw If the argument key value isn't found, then the function throws an error
         */
        hur_nodiscard inline const Type &at(const keyType &keyval) const;

        /**
         * \brief Finds or inserts an element with the specified key.
         * \param[in] keyval - The key value to find or insert.
         * \return A reference to the data value of the inserted element.
         * \remarks	If the argument key value isn't found, then it's inserted along with the default value of the data type.
         */
        inline Type &operator[](const keyType &keyval);

        /**
         * \brief Finds or inserts an element with the specified key.
         * \param[in] keyval - The key value to find or insert.
         * \return A reference to the data value of the inserted element.
         * \remarks	If the argument key value isn't found, then it's inserted along with the default value of the data type.
         */
        inline Type &operator[](keyType &&keyval);

        /**
         * \brief inds an element in a unordered_map with a specified key value.
         * \param[in] keyval - The key value to find.
         * \return A reference to the data value of the element found.
         * \throw If the argument key value isn't found, then the function throws an error
         */
        inline const Type &operator[](const keyType &keyval) const;
    };
} // namespace OpenHurricane
#include "HashMap.inl"