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

SparseDataMatrix::SparseDataMatrix() : nrows{0}, ncols{0}, data{}, colIndex{}, rowPtr{0} {}

SparseDataMatrix::SparseDataMatrix(size_t nrows, size_t ncols)
    : nrows{nrows}, ncols{ncols}, data{}, colIndex{}, rowPtr(ncols, 0) {}

size_t SparseDataMatrix::getNrows() const { return nrows; }

size_t SparseDataMatrix::getNcols() const { return ncols; }

void SparseDataMatrix::resize(size_t nrows, size_t ncols) {
  const auto oldNrows = this->nrows;
  const auto oldNcols = this->ncols;

  const auto newNrows = nrows;
  const auto newNcols = ncols;

  /**
   *  rows
   */

  // add new empty rows
  if (oldNrows < newNrows) {
    const auto lastRow = rowPtr.back();
    rowPtr.resize(newNrows, lastRow);

    // eliminate last rows
  } else {
    // the last elements:
    const auto lastElem = rowPtr[newNrows];
    rowPtr.resize(newNrows);
    colIndex.resize(lastElem);
    data.resize(lastElem);
  }

  /**-
   * columns
   */

  // resizing only required, if new size is smaller then old one.
  if (oldNcols > newNcols) {
    // for each row
    for (auto i = 0u; i < newNrows; i++) {
      // determine last element in current row
      const auto upper = i + 1 < rowPtr.size() ? rowPtr[i + 1] : colIndex.size();
      // for each element in current row
      for (auto j = rowPtr[i]; j < upper; j++) {
        if (colIndex[j] >= newNcols) {
          const auto shift = upper - j - 1;  // index!

          // erase the actual values;
          data.erase(data.begin() + j, data.begin() + upper);
          // erase the column indices
          colIndex.erase(colIndex.begin() + j, colIndex.begin() + upper);
          // shift the row offset
          for (auto idx = rowPtr.begin() + i; idx < rowPtr.end(); idx++) {
            *idx -= shift;
          }
        }
      }
    }
  }

  // finally adapt matrix size:
  this->nrows = newNrows;
  this->ncols = newNcols;
}

const std::vector<size_t>& SparseDataMatrix::getColIndexVector() const { return colIndex; }

const std::vector<double>& SparseDataMatrix::getDataVector() const { return data; }

const std::vector<size_t>& SparseDataMatrix::getRowPtrVector() const { return rowPtr; }

std::vector<size_t>& SparseDataMatrix::getColIndexVector() {
  return const_cast<std::vector<size_t>&>(
      static_cast<const SparseDataMatrix&>(*this).getColIndexVector());
}

std::vector<double>& SparseDataMatrix::getDataVector() {
  return const_cast<std::vector<double>&>(
      static_cast<const SparseDataMatrix&>(*this).getDataVector());
}

std::vector<size_t>& SparseDataMatrix::getRowPtrVector() {
  return const_cast<std::vector<size_t>&>(
      static_cast<const SparseDataMatrix&>(*this).getRowPtrVector());
}

void SparseDataMatrix::fromDataMatrix(const DataMatrix& in, SparseDataMatrix& out,
                                      double threshold) {
  auto inRows = in.getNrows();
  auto inCols = in.getNcols();

  out.resize(inRows, inCols);

  auto tmpIn = 0.0;
  for (auto i = 0u; i < inRows; i++) {
    out.getRowPtrVector()[i] = out.getDataVector().size();
    for (auto j = 0u; j < inCols; j++) {
      tmpIn = in.get(i, j);
      if (std::abs(tmpIn) > threshold) {
        out.getDataVector().push_back(tmpIn);
        out.getColIndexVector().push_back(j);
      }
    }
  }
}

void SparseDataMatrix::toDataMatrix(const SparseDataMatrix& in, DataMatrix& out) {
  auto inRows = in.getNrows();
  auto inCols = in.getNcols();
  const auto& inRowPtr = in.getRowPtrVector();
  const auto& inIdx = in.getColIndexVector();

  out.resize(inRows, inCols);
  out.setAll(0.0);

#pragma omp parallel for
  for (auto i = 0u; i < inRows; i++) {
    const auto upper = i + 1 < inRowPtr.size() ? inRowPtr[i + 1] : inIdx.size();
    for (auto j = inRowPtr[i]; j < upper; j++) {
      out.set(i, inIdx[j], in.getDataVector()[j]);
    }
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
