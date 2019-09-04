// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * DBMatOffline specialization as a base class for all algorithms based on gaussian elimination on
 * a dense matrix.
 */
class DBMatOfflineGE : public DBMatOffline {
 public:
  explicit DBMatOfflineGE(const std::string& fileName);

  /**
   * Builds the right hand side matrix with identity regularization term
   * @param grid The grid object the matrix is based on
   * @param regularizationConfig Configures the regularization which is incorporated into the lhs
   */
  void buildMatrix(Grid* grid, RegularizationConfiguration& regularizationConfig) override;

 protected:
  DBMatOfflineGE();
};

} /* namespace datadriven */
} /* namespace sgpp */
