/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseDataMatrix.cpp
 *
 *  Created on: Feb 4, 2017
 *      Author: Michael Lettrich
 */

#include "SparseDataMatrix.hpp"

#include <sgpp/base/exception/application_exception.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

SparseDataMatrix::SparseDataMatrix() : nrows{0}, ncols{0}, data{}, colIndex{}, rowIndicator{} {}

SparseDataMatrix::SparseDataMatrix(size_t nrows, size_t ncols)
    : nrows{nrows}, ncols{ncols}, data{}, colIndex{}, rowIndicator(ncols) {}

size_t SparseDataMatrix::getNrows() const { return nrows; }

size_t SparseDataMatrix::getNcols() const { return ncols; }

void SparseDataMatrix::resize(size_t nrows, size_t ncols) {
  const std::string msg{__PRETTY_FUNCTION__, " : not implemented yet.\n"};
  throw(msg.c_str());
}

const std::vector<size_t>& SparseDataMatrix::getColIndexVector() const { return colIndex; }

const std::vector<double>& SparseDataMatrix::getDataVector() const { return data; }

const std::vector<size_t>& SparseDataMatrix::getRowIndicatorVector() const { return rowIndicator; }

std::vector<size_t>& SparseDataMatrix::getColIndexVector() {
  return const_cast<std::vector<size_t>&>(
      static_cast<const SparseDataMatrix&>(*this).getColIndexVector());
}

std::vector<double>& SparseDataMatrix::getDataVector() {
  return const_cast<std::vector<double>&>(
      static_cast<const SparseDataMatrix&>(*this).getDataVector());
}

std::vector<size_t>& SparseDataMatrix::getRowIndicatorVector() {
  return const_cast<std::vector<size_t>&>(
      static_cast<const SparseDataMatrix&>(*this).getRowIndicatorVector());
}

void SparseDataMatrix::fromDataMatrix(const DataMatrix& in, SparseDataMatrix& out,
                                      double threshold) {
  auto rows = in.getNrows();
  auto cols = in.getNcols();

  out.resize(rows, cols);

  auto tmpIn = 0.0;
  // #pragma omp parallel for private(tmpIn);
  for (auto i = 0u; i < rows; i++) {
    out.getRowIndicatorVector()[i] = out.getDataVector().size();
    for (auto j = 0u; j < cols; j++) {
      tmpIn = in.get(i, j);
      if (tmpIn > threshold) {
        out.getDataVector().push_back(tmpIn);
        out.getColIndexVector().push_back(j);
      }
    }
  }
}

void SparseDataMatrix::toDataMatrix(const SparseDataMatrix& in, DataMatrix& out) {
  auto rows = in.getNrows();
  auto cols = in.getNcols();

  const auto& rowVec = in.getRowIndicatorVector();

  out.resize(rows, cols);
  out.setAll(0.0);

  // #pragma omp parallel for private(tmpIn);
  for (auto i = 0u; i < in.getNrows(); i++) {
    const auto upper = i + 1 > rowVec.size() ? rowVec[i + 1] : in.getColIndexVector().size();
    for (auto j = rowVec[i]; j < upper; j++) {
      out.set(i, j, in.getDataVector()[i]);
    }
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
