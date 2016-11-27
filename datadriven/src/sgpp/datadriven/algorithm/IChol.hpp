/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * IChol.hpp
 *
 *  Created on: Nov 26, 2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;

class IChol {
 public:
  /**
   * decompose in place
   */
  static void decompose(DataMatrix& matrix, size_t sweeps);
  /**
   * decomposition is stored in results matrix
   */
  static void decompose(const DataMatrix& matrix, DataMatrix& result, size_t sweeps);
  /**
   * do a cholesky update on the last n rows
   */
  static void updateLastNRows(DataMatrix& matrix, size_t numRows, size_t sweeps);
};

} /* namespace datadriven */
} /* namespace sgpp */
