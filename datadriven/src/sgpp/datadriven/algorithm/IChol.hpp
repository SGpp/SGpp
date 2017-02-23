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

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/algorithm/SparseDataMatrix.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataVector;

class IChol {
 public:
  /**
   * decompose in place
   */
  static void decompose(SparseDataMatrix& matrix, size_t sweeps);

  /**
   * do a cholesky update on the last n rows
   */
  static void updateLastNRows(SparseDataMatrix& matrix, size_t numRows, size_t sweeps);

  static void normToUnitDiagonal(SparseDataMatrix& matrix, DataVector& norms);

  static void reaplyDiagonal(SparseDataMatrix& matrix, DataVector& norms);
};

} /* namespace datadriven */
} /* namespace sgpp */
