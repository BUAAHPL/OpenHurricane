/*!
 * \file KDTree.hpp
 * \brief Header of binaryTree.
 *       The subroutines and functions are in the <i>KDTree.cpp</i> file.
 * \author Yang Hongzhen
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

#include "Lists.hpp"
#include "dataStructure.hpp"

namespace OpenHurricane {

    
    template <class Node, class Type> class binaryTree {
    protected:
        Type *root_;

        integer *index_;

        integer *splitDim_;

        Node *parent_;

        Node *leftChild_;

        Node *rightChild_;

    public:
        inline binaryTree()
            : root_(nullptr), index_(nullptr), splitDim_(nullptr), parent_(nullptr),
              leftChild_(nullptr), rightChild_(nullptr) {}

        inline virtual ~binaryTree() noexcept { clear(this); }

        hur_nodiscard inline bool isEmpty() const noexcept { return root_ == nullptr; }

        hur_nodiscard inline bool isRoot() const noexcept {
            return (root_ != nullptr) && parent_ == nullptr;
        }

        hur_nodiscard inline bool isLeft() const noexcept {
            return parent_->leftChild_->root_ == root_;
        }

        hur_nodiscard inline bool isRight() const noexcept {
            return parent_->rightchild_->root_ == root_;
        }

        hur_nodiscard inline bool isLeaf() const noexcept {
            return (root_ != nullptr) && rightChild_ == nullptr && leftChild_ == nullptr;
        }

        inline void setRoot(const Type root) noexcept { root_ = new Type(root); }

        inline void setParent(const Node *parent) noexcept { parent_ = parent; }

        inline void setLeft(const Node *left) noexcept { leftChild_ = left; }

        inline void setRight(const Node *right) noexcept { rightChild_ = right; }

        void clear(binaryTree *tree) noexcept {
            if (tree == nullptr) {
                return;
            }
            if (tree->leftChild_ != nullptr) {
                clear(tree->leftChild_);
            }

            if (tree->rightChild_ != nullptr) {
                clear(tree->rightChild_);
            }

            if (tree != nullptr) {
                HurDelete(tree->root_);
                HurDelete(tree->index_);
                HurDelete(tree->splitDim_);
                HurDelete(tree->leftChild_);
                HurDelete(tree->rightChild_);
            }
            return;
        }
    };
    template <class Node> using vectorTree = binaryTree<Node, vector>;

    class KDTree : public vectorTree<KDTree> {
    private:
    public:
        KDTree();

        inline ~KDTree() noexcept {};

        hur_nodiscard real middleValue(const integerList &nodeset, const vectorList &nodes,
                                       integer depth) const;

        void buildKDTree(KDTree *tree, const integerList &nodeset, const vectorList &nodes,
                         integer depth);

        integerList getNearestNbr(integer id, const vector goal, KDTree *tree);

        void printTree(KDTree *tree, integer depth);

        void checkTree(KDTree *tree, integer depth);
    };
} // namespace OpenHurricane