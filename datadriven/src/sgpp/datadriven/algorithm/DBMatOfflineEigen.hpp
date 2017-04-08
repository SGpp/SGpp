/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineEigen.hpp
 *
 *  Created on: 02.03.2017
 *      Author: Michael Lettrich
 */
#ifdef USE_GSL

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <string>

namespace sgpp {
namespace datadriven {

class DBMatOfflineEigen : public DBMatOffline {
 public:
  /**
   * Constructor
   *
   * @param oc configuration for this offline object
   */
  explicit DBMatOfflineEigen(const DBMatDensityConfiguration& oc);

  explicit DBMatOfflineEigen(const std::string& fileName);

  DBMatOffline* clone() override;

  /**
   * Decomposes the matrix according to the chosen decomposition type.
   * The number of rows of the stored result depends on the decomposition type.
   */
  void decomposeMatrix() override;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* USE_GSL */
