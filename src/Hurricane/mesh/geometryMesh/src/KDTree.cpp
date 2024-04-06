/*!
 * \file KDTree.inl
 * \brief The subroutines and functions of CFD time advance iteration
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

#include "KDTree.hpp"
#include "stack"

OpenHurricane::KDTree::KDTree() : vectorTree<KDTree>() {}

hur_nodiscard OpenHurricane::real OpenHurricane::KDTree::middleValue(
    const integerList &nodeset, const vectorList &nodes, integer depth) const {
    integer k = feature<vector>::nElements_;
    integer splitAttribute = depth % k;
    realList dimCor(nodeset.size());

    for (integer i = 0; i < nodeset.size(); i++) {
        integer index = nodeset[i];
        dimCor[i] = nodes[index][splitAttribute];
    }
    std::sort(dimCor.data(), dimCor.data() + dimCor.size());
    integer pos = nodeset.size() / 2;
    return real(dimCor[pos]);
}

void OpenHurricane::KDTree::buildKDTree(KDTree *tree,
                                        const integerList &nodeset,
                                        const vectorList &nodes,
                                        integer depth) {
    if (nodeset.size() == 0) {
        return;
    }
    if (nodeset.size() == 1) {
        tree->root_ = new vector(nodes[nodeset[0]]);
        tree->index_ = new integer(nodeset[0]);
        return;
    }

    integer k = feature<vector>::nElements_;

    integer splitAttribute = depth % k;

    double splitValue = middleValue(nodeset, nodes, depth);

    integer left = 0;
    integer right = 0;

    for (integer i = 0; i < nodeset.size(); ++i) {
        integer index = nodeset[i];
        if (nodes[index][splitAttribute] == splitValue) {
            if (tree->root_ == nullptr) {
                tree->root_ = new vector(nodes[index]);
                tree->index_ = new integer(index);
                tree->splitDim_ = new integer(splitAttribute);
            } else {
                left++;
            }
        } else {
            if (nodes[index][splitAttribute] < splitValue) {
                left++;
            } else {
                right++;
            }
        }
    }

    integerList leftSet(left);
    integerList rightSet(right);
    left = 0;
    right = 0;

    for (integer i = 0; i < nodeset.size(); ++i) {
        integer index = nodeset[i];
        if (index != *tree->index_)
        {
            if (nodes[index][splitAttribute] <= splitValue) {
                leftSet[left++] = index;
            } else 
            {
                rightSet[right++] = index;
            }
        }
    }

    if (leftSet.size() != 0) {
        tree->leftChild_ = new KDTree;
        tree->leftChild_->parent_ = tree;
        buildKDTree(tree->leftChild_, leftSet, nodes, depth + 1);
    }

    if (rightSet.size() != 0) {
        tree->rightChild_ = new KDTree;
        tree->rightChild_->parent_ = tree;
        buildKDTree(tree->rightChild_, rightSet, nodes, depth + 1);
    }
}

OpenHurricane::integerList
OpenHurricane::KDTree::getNearestNbr(integer id, const vector goal,
                                     KDTree *tree) {

    KDTree *currentTree = tree;
    vector currentNearest = *currentTree->root_;
    integer currentIndex = *currentTree->index_;
    std::stack<KDTree *> searchPath;
    while (!currentTree->isLeaf())
    {
        searchPath.push(currentTree);
        integer index = *currentTree->splitDim_;
        if (currentTree->rightChild_ == nullptr ||
            goal[index] <= (*currentTree->root_)[index])
        {
            currentTree = currentTree->leftChild_;
        } else {
            currentTree = currentTree->rightChild_;
        }
    }

    currentNearest = *currentTree->root_;
    currentIndex = *currentTree->index_;
    real currentDistance = (currentNearest - goal).magnitude();

    while (searchPath.size() != 0) {
        KDTree *searchDistrict = searchPath.top();
        searchPath.pop();

        if (searchDistrict->isLeaf()) {
            real distance = (*searchDistrict->root_ - goal).magnitude();
            if (distance < currentDistance) {
                currentDistance = distance;
                currentTree = searchDistrict;
                currentNearest = *currentTree->root_;
                currentIndex = *currentTree->index_;
            }
        } else {
            integer index = *searchDistrict->splitDim_;
            real districtDistance =
                fabs(goal[index] - (*searchDistrict->root_)[index]);

            if (districtDistance < currentDistance) {
                real distance = (goal - *searchDistrict->root_).magnitude();

                if (distance < currentDistance) {
                    currentDistance = distance;
                    currentTree = searchDistrict;
                    currentNearest = *currentTree->root_;
                    currentIndex = *currentTree->index_;
                }

                KDTree *addSearch(nullptr);
                if (goal[index] <= (*searchDistrict->root_)[index]) {
                    addSearch = searchDistrict->rightChild_;
                } else {
                    addSearch = searchDistrict->leftChild_;
                }

                if (addSearch != nullptr) {
                    searchPath.push(addSearch);
                }
            }
        }
    }
    integerList nearest(1, currentIndex);
    return nearest;
}

void OpenHurricane::KDTree::printTree(KDTree *tree, integer depth) {
    if (!tree->isLeaf()) {
        std::cout << "depth = " << depth << " root = " << *tree->index_;

        if (tree->leftChild_ != nullptr) {
            std::cout << "  left = " << *tree->leftChild_->index_;
        }
        if (tree->rightChild_ != nullptr) {
            std::cout << "  right = " << *tree->rightChild_->index_;
        }
        std::cout << "" << std::endl;
        if (tree->leftChild_ != nullptr) {
            printTree(tree->leftChild_, depth + 1);
        }
        if (tree->rightChild_ != nullptr) {
            printTree(tree->rightChild_, depth + 1);
        }
    }
}

void OpenHurricane::KDTree::checkTree(KDTree *tree, integer depth) {
    if (!tree->isLeaf()) {
        integer dim = depth % 3;
        if (*tree->splitDim_ != dim) {
        }
        if ((*tree->leftChild_->root_)[dim] > (*tree->root_)[dim]) {
        }
        if (tree->rightChild_ != nullptr) {
            if ((*tree->rightChild_->root_)[dim] <= (*tree->root_)[dim]) {
            }
        }
        if (tree->parent_ != nullptr) {
            if (tree->isLeft()) {
                integer id = *tree->parent_->splitDim_;
                if (tree->rightChild_ != nullptr) {
                    if ((*tree->rightChild_->root_)[id] >
                        (*tree->parent_->root_)[id]) {
                    }
                } else {
                    if ((*tree->leftChild_->root_)[id] >
                        (*tree->parent_->root_)[id]) {
                    }
                }
            } else {
                integer id = *tree->parent_->splitDim_;
                if ((*tree->leftChild_->root_)[id] <=
                    (*tree->parent_->root_)[id]) {
                }
            }
        }

        checkTree(tree->leftChild_, depth + 1);

        if (tree->rightChild_ != nullptr) {
            checkTree(tree->rightChild_, depth + 1);
        }
    }
}
