/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineDenseIChol.hpp
 *
 *  Created on: Apr 15, 2017
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>

namespace sgpp {
namespace datadriven {

class DBMatOfflineDenseIChol : public DBMatOfflineChol {
 public:
  /**
   * Constructor
   *
   * @param oc configuration for this offline object
   */
  explicit DBMatOfflineDenseIChol(const DBMatDensityConfiguration& oc);

  explicit DBMatOfflineDenseIChol(const std::string& fileName);

  DBMatOffline* clone() override;

  /**
   * Decomposes the matrix according to the chosen decomposition type.
   * The number of rows of the stored result depends on the decomposition type.
   */
  void decomposeMatrix() override;

  /**
   * Updates offline cholesky factorization based on coarsed (deletedPoints)
   * and refined (newPoints) gridPoints
   *
   * @param deletedPoints list of indices of last coarsed points
   * @param newPoints amount of refined points
   */
  void choleskyModification(size_t newPoints, std::list<size_t> deletedPoints,
                            double lambda) override;

 protected:
  void ichol(const DataMatrix& matrix, DataMatrix& result, size_t sweeps = 5, size_t startRow = 0);
};

} /* namespace datadriven */
} /* namespace sgpp */
