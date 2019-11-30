// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/tools/vpTree/VpNode.hpp>
#include <sgpp/datadriven/tools/vpTree/VpHeapItem.hpp>

#include <vector>
#include<utility>
#include <queue>
namespace sgpp {
namespace datadriven {
/**
 * Class for the Vantage Point Tree.
 * Based on the code by Steve Hanov's great tutorial
 * at http://stevehanov.ca/blog/index.php?id=130
 */
class VpTree {
 public:
  /**
   * Constructor
   * @param matrix
   * Matrix of points to use
   */
  explicit VpTree(sgpp::base::DataMatrix matrix);

  // Destructor
  ~VpTree() {
    delete root;
  }

  std::priority_queue<VpHeapItem>  getNearestNeighbors(sgpp::base::DataVector &target,
      size_t noNearestNeighbors);

  void update(sgpp::base::DataMatrix &matrix);

  void updateHard(sgpp::base::DataMatrix &matrix);

  sgpp::base::DataMatrix &getStoredItems();

  /**
  * Defines the euclidean distance metric between two points
  */
  static double euclideanDistance(sgpp::base::DataVector point1, sgpp::base::DataVector point2);

  void printPreorder();

 private:
  /**
   * Node that defines the root of the tree
   */
  VpNode* root;
  /**
   * Value that keeps track the distance of the latest found nearest neighbor
   */
  double tau;

  /**
   * Matrix to keep track of all stored points when building the tree
   */
  sgpp::base::DataMatrix storedItems;

  /**
   * Method that builds a VpTree recursevily
   * @param startIndex Index of the matrix of the points to store
   * @param endIndex End index of the matrix of the poinst to store
   * @return The node corresponding to the roort of the tree
   */
  VpNode* buildRecursively(size_t startIndex, size_t endIndex);

  void insertNewNode(size_t indexNewPoint);

  VpNode* findInsertionNode(VpNode* root, sgpp::base::DataVector &newPoint);

  void searchRecursively(VpNode* &node, sgpp::base::DataVector &target,
      size_t noNearestNeighbors, std::priority_queue<VpHeapItem> &heap);


  void swap(size_t index1, size_t index2);

  std::vector<std::pair <size_t, double>>  getDistances(size_t startIndex, size_t endIndex);

  void sortByDistances(size_t index1, size_t index2);

  void printPreorder(VpNode* node);

};
}  // namespace datadriven
}  // namespace sgpp
