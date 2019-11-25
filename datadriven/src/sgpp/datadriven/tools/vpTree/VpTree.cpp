// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/vpTree/VpTree.hpp>
#include <cfloat>
#include <utility>
#include <algorithm>
#include <vector>
#include <queue>
#include <iostream>

using sgpp::base::DataVector;
using sgpp::base::DataMatrix;
using sgpp::datadriven::VpNode;

namespace sgpp {
namespace datadriven {

VpTree::VpTree(DataMatrix matrix) {
  this->tau = DBL_MAX;
  this->storedItems = matrix;
  this->root = buildRecursively(0, this->storedItems.getNrows());
}


VpTree::~VpTree() {
  delete root;
}

// Function that uses the tree to find the k nearest neighbors of target
std::priority_queue<VpHeapItem> VpTree::getNearestNeighbors(DataVector &target,
    size_t noNearestNeighbors) {
  std::priority_queue<VpHeapItem> heap;

  if (noNearestNeighbors >= storedItems.getNrows()) {
    noNearestNeighbors = storedItems.getNrows();
  }
  tau = DBL_MAX;

  searchRecursively(root, target, noNearestNeighbors, heap);

  return heap;
}

void VpTree::searchRecursively(VpNode* &node, DataVector &target, size_t noNearestNeighbors,
    std::priority_queue<VpHeapItem> &heap) {
  if (node == nullptr) {
    return;
  }

  DataVector currentVector(target.size());

  storedItems.getRow(node->index, currentVector);

  double distance = euclideanDistance(currentVector, target);

  if (distance < tau && distance != 0) { // the second condition is to skip the same point
    if (heap.size() == noNearestNeighbors) {
      heap.pop();
    }
    heap.push(VpHeapItem(node->index, distance));
    if (heap.size() == noNearestNeighbors) {
      tau = heap.top().distance;
    }
  }

  // Return if we arrived at a leaf
  if (node->left == NULL && node->right == NULL) {
    return;
  }

  // If the target lies within the radius of ball
  if (distance <= node->threshold) {
    if (distance - tau <= node->threshold) {
      // if there can still be neighbors inside the ball, recursively search left child first
      searchRecursively(node->left, target, noNearestNeighbors, heap);
    }

    if (distance + tau >= node->threshold) {
      // if there can still be neighbors outside the ball, recursively search right child
      searchRecursively(node->right, target, noNearestNeighbors, heap);
    }
    // If the target lies outsize the radius of the ball
  } else {
    if (distance + tau >= node->threshold) {
      // if there can still be neighbors outside the ball, recursively search right child first
      searchRecursively(node->right, target, noNearestNeighbors, heap);
    }

    if (distance - tau <= node->threshold) {
      // if there can still be neighbors inside the ball, recursively search left child
      searchRecursively(node->left, target, noNearestNeighbors, heap);
    }
  }
}

void VpTree::update(DataMatrix &matrix) {
  DataVector newRow(matrix.getNcols());
  size_t lastIndex = storedItems.getNrows();
  for (size_t index = 0; index < matrix.getNrows() ; index++) {
    matrix.getRow(index, newRow);
    storedItems.appendRow(newRow);
    insertNewNode(lastIndex + index);
  }

}

void VpTree::updateHard(sgpp::base::DataMatrix &matrix) {
  DataVector newRow(matrix.getNcols());
  size_t lastIndex = storedItems.getNrows();
  for (size_t index = 0; index < matrix.getNrows() ; index++) {
    matrix.getRow(index, newRow);
    storedItems.appendRow(newRow);
  }
  root = buildRecursively(0, storedItems.getNrows());
}
VpNode* VpTree::buildRecursively(size_t startIndex, size_t endIndex)  {
    if (endIndex == startIndex) {     // indicates that we're done here!
      return nullptr;
    }
    // startIndex index is center of current VpNode
    DataVector vector1(storedItems.getNcols());
    storedItems.getRow(startIndex, vector1);
    VpNode* node = new VpNode();
    node->index = startIndex;
    node->threshold = 0.0;


    if (endIndex - startIndex > 1) {
      size_t median = (endIndex + startIndex) / 2;

      sortByDistances(startIndex, endIndex);

      DataVector vector2(storedItems.getNcols());
      storedItems.getRow(median, vector2);

      node->threshold = euclideanDistance(vector1, vector2);

      node->left = buildRecursively(startIndex + 1, median+1);
      node->right = buildRecursively(median+1, endIndex);
    }

    return node;
}

void VpTree::insertNewNode(size_t indexNewPoint) {
  DataVector newPoint(storedItems.getNcols());
  storedItems.getRow(indexNewPoint, newPoint);

  VpNode* newNode = new VpNode();
  newNode->index = indexNewPoint;
  newNode->threshold = 0;

  VpNode* insertionNode = findInsertionNode(root, newPoint);

  DataVector vantagePoint(storedItems.getNcols());
  storedItems.getRow(insertionNode->index, vantagePoint);

  double distance = euclideanDistance(newPoint, vantagePoint);

  if (insertionNode->threshold == 0) {
    insertionNode->threshold = distance;
    insertionNode->left = newNode;
  } else {
    if (distance <= insertionNode->threshold) {
      insertionNode->left = newNode;
    } else {
      insertionNode->right = newNode;
    }
  }
}

VpNode* VpTree::findInsertionNode(VpNode* node, DataVector &newPoint) {
  if (node->threshold == 0) {
    return node;
  }
  DataVector vantagePoint(storedItems.getNcols());

  storedItems.getRow(node->index, vantagePoint);

  double distance = euclideanDistance(newPoint, vantagePoint);

  if (distance <= node->threshold) {
    if (node->left == nullptr) {
      return node;
    } else {
      return findInsertionNode(node->left, newPoint);
    }
  } else {
    if (node->right == nullptr) {
      return node;
    } else {
      return findInsertionNode(node->right, newPoint);
    }
  }
}

double VpTree::euclideanDistance(DataVector point1, DataVector point2) {
  point1.sub(point2);
  return point1.l2Norm();
}

void VpTree::swap(size_t index1, size_t index2) {
  DataVector vector1(storedItems.getNcols());
  storedItems.getRow(index1, vector1);

  DataVector vector2(storedItems.getNcols());
  storedItems.getRow(index2, vector2);

  storedItems.setRow(index1, vector2);
  storedItems.setRow(index2, vector1);
}

std::vector<std::pair <size_t, double>>  VpTree::getDistances(size_t startIndex, size_t endIndex) {
  DataVector vantagePoint(storedItems.getNcols());
  storedItems.getRow(startIndex, vantagePoint);
  DataVector currentRow(storedItems.getNcols());

  std::vector<std::pair <size_t, double>> distances;

  for (size_t index = startIndex+1; index < endIndex ; index++) {
    storedItems.getRow(index, currentRow);
    double distance = euclideanDistance(vantagePoint, currentRow);
    distances.push_back(std::make_pair(index, distance));
  }
  return distances;
}

void VpTree::sortByDistances(size_t index1, size_t index2) {
  auto distances = getDistances(index1, index2);

  std::sort(distances.begin(), distances.end(), [](const std::pair<size_t, double> &pair1,
      const std::pair<size_t, double> &pair2) -> bool {
    return pair1.second < pair2.second;
  });

  std::vector<size_t> visitedIndexes;

  for (size_t index = 0; index< index2-index1-1; index++) {
    if (visitedIndexes.size()==index2-index1-1) {
      return;
    }
    if (std::find(visitedIndexes.begin(), visitedIndexes.end(), distances[index].first) != visitedIndexes.end()) {
      continue;
    }
    visitedIndexes.push_back(index1+index+1);
    visitedIndexes.push_back(distances[index].first);
    swap(index1+index+1, distances[index].first);
  }
}

DataMatrix &VpTree::getStoredItems() {
  return this->storedItems;
}

void VpTree::printPreorder() {
  printPreorder(root);
}
void VpTree::printPreorder(VpNode* node) {

  if(node == nullptr) {
    return;
  }
  DataVector currentNode(storedItems.getNcols());
  storedItems.getRow(node->index, currentNode);
  std::cout << currentNode.toString() << ", "
  << std::to_string(node->threshold) << std::endl;

  printPreorder(node->left);
  printPreorder(node->right);

}
}  // namespace datadriven
}  // namespace sgpp
