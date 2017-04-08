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

using sgpp::base::DataVector;
using sgpp::base::DataMatrix;
using sgpp::base::Grid;

class DBMatOfflineIChol : public DBMatOfflineChol {
 public:
  /**
   * Constructor
   *
   * @param oc configuration for this offline object
   */
  explicit DBMatOfflineIChol(const DBMatDensityConfiguration& oc);

  explicit DBMatOfflineIChol(const std::string& fileName);

  DBMatOffline* clone() override;

  /**
   * Decomposes the matrix according to the chosen decomposition type.
   * The number of rows of the stored result depends on the decomposition type.
   */
  virtual void decomposeMatrix() override;
};

} /* namespace datadriven */
} /* namespace sgpp */
