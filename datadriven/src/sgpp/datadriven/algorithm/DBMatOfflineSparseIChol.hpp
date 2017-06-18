/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineIChol.hpp
 *
 *  Created on: Feb 27, 2017
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>

namespace sgpp {
namespace datadriven {

/**
 * TODO(lettrich) : write documentation
 */

class DBMatOfflineSparseIChol : public DBMatOfflineChol {
 public:
  /**
   * Constructor
   *
   * @param oc configuration for this offline object
   */
  explicit DBMatOfflineSparseIChol(const DBMatDensityConfiguration& oc);

  explicit DBMatOfflineSparseIChol(const std::string& fileName);

  DBMatOffline* clone() override;

  /**
   * Decomposes the matrix according to the chosen decomposition type.
   * The number of rows of the stored result depends on the decomposition type.
   */
  void decomposeMatrix() override;
};

} /* namespace datadriven */
} /* namespace sgpp */
