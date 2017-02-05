/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseDataMatrix.hpp
 *
 *  Created on: Feb 4, 2017
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;

class SparseDataMatrix {
 public:
  SparseDataMatrix();
  SparseDataMatrix(size_t nrows, size_t ncols);

  size_t getNrows() const;

  size_t getNcols() const;

  void resize(size_t nrows, size_t ncols);

  const std::vector<size_t>& getColIndexVector() const;

  const std::vector<double>& getDataVector() const;

  const std::vector<size_t>& getRowPtrVector() const;

  std::vector<size_t>& getColIndexVector();

  std::vector<double>& getDataVector();

  std::vector<size_t>& getRowPtrVector();

  static void fromDataMatrix(const DataMatrix& in, SparseDataMatrix& out,
                             double threshold = std::numeric_limits<double>::epsilon());

  static void toDataMatrix(const SparseDataMatrix& in, DataMatrix& out);

 private:
  size_t nrows;
  size_t ncols;

  std::vector<double> data;
  std::vector<size_t> colIndex;
  std::vector<size_t> rowPtr;
};

} /* namespace datadriven */
} /* namespace sgpp */
