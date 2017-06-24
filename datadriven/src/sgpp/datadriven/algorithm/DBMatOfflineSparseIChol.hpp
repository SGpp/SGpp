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

#include <sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * DBMatOfflineChol specialization that uses a parallel, iterative incomplete cholesky factorization
 * on a a sparse matrix. The current implementation is a proof of concept.
 */

class DBMatOfflineSparseIChol : public DBMatOfflineDenseIChol {
 public:
  explicit DBMatOfflineSparseIChol(const DBMatDensityConfiguration& oc);

  explicit DBMatOfflineSparseIChol(const std::string& fileName);

  DBMatOffline* clone() override;

  void decomposeMatrix() override;

  /**
   * the actual incomplete Cholesky factorization using sparse matrices internally. This
   * implementation is an experimental POC and should not be used in production.
   */
  static void ichol(const DataMatrix& matrix, DataMatrix& result, size_t sweeps = 4,
                    size_t startRow = 0);
};

} /* namespace datadriven */
} /* namespace sgpp */
