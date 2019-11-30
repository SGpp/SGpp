// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace sgpp {
namespace datadriven {
/**
 * Class that defines the nodes of the Vantage Point Tree based on the DataVector structure
 * of the SG++ Toolbox
 */
struct VpNode {
  VpNode() = default;

  VpNode(VpNode &&rhs) = default;

  VpNode(const VpNode &rhs) = default;

  VpNode &operator=(VpNode &&rhs) = default;
  /**
   * Destructor.
   */
  ~VpNode() {
    delete left;
    delete right;
  }

  /* Pointer to the node with points closer than the threshold */
  VpNode* left;

  /* Pointer to the node with points farther than the threshold */
  VpNode* right;

  /* Index of the point of used to keep track of it */
  size_t index;
  /* Threshold that divides the tree into 2 */
  double threshold;

};
}  // namespace datadriven
}  // namespace sgpp

