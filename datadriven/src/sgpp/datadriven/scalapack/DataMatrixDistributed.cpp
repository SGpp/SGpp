/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataMatrixDistributed.cpp
 *
 * Created on: Jan 14, 2019
 *     Author: Jan Schopohl
 */
#ifdef USE_SCALAPACK

#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>

#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>
#include <sgpp/datadriven/scalapack/blacs.hpp>

#include <iostream>
#include <sgpp/base/exception/not_implemented_exception.hpp>

namespace sgpp {
namespace datadriven {

DataMatrixDistributed::DataMatrixDistributed(std::shared_ptr<BlacsProcessGrid> grid, int globalRows,
                                             int globalColumns, DataMatrixDistributed::DTYPE dtype,
                                             int columnBlockSize, int rowBlockSize, double value)
    : grid(grid),
      globalRows(globalRows),
      globalColumns(globalColumns),
      rowBlockSize(rowBlockSize),
      columnBlockSize(columnBlockSize) {
  descriptor[dtype_] = static_cast<int>(dtype);

  // TODO(jan) implement other matrix types
  if (dtype != DataMatrixDistributed::DTYPE::DENSE) {
    throw sgpp::base::not_implemented_exception("Only Dense Matrix implemented at the moment");
  }

  localRows = numroc_(globalRows, rowBlockSize, grid->getCurrentRow(), 0, grid->getTotalRows());
  localColumns =
      numroc_(globalColumns, columnBlockSize, grid->getCurrentColumn(), 0, grid->getTotalColumns());

  // initialize local matrix
  localData.assign(localRows * localColumns, value);

  descriptor[ctxt_] = grid->getContextHandle();
  descriptor[m_] = globalRows;
  descriptor[n_] = globalColumns;
  descriptor[mb_] = rowBlockSize;
  descriptor[nb_] = columnBlockSize;
  descriptor[mb_] = localRows;
  descriptor[nb_] = localColumns;
  descriptor[rsrc_] = 0;
  descriptor[csrc_] = 0;
  descriptor[lld_] = localRows;
}

DataMatrixDistributed::DataMatrixDistributed(double* input, std::shared_ptr<BlacsProcessGrid> grid,
                                             int globalRows, int globalColumns,
                                             DataMatrixDistributed::DTYPE dtype,
                                             int columnBlockSize, int rowBlockSize)
    : DataMatrixDistributed(grid, globalRows, globalColumns, dtype, columnBlockSize, rowBlockSize,
                            0.0) {
  int desca[dlen_];
  desca[ctxt_] = -1;
  BlacsProcessGrid localGrid = BlacsProcessGrid(1, 1);
  if (grid->getCurrentProcess() == 0) {
    desca[ctxt_] = localGrid.getContextHandle();
    desca[m_] = globalRows;
    desca[n_] = globalColumns;
    desca[mb_] = globalRows;
    desca[nb_] = globalColumns;
    desca[rsrc_] = 0;
    desca[csrc_] = 0;
    desca[lld_] = globalRows;
  }
  pdgemr2d_(globalRows, globalColumns, input, 1, 1, desca, this->getLocalPointer(), 1, 1,
            this->descriptor, this->grid->getContextHandle());
}

DataMatrix DataMatrixDistributed::toLocalDataMatrix() const {
  DataMatrix localMatrix(globalRows, globalColumns);

  int desca[dlen_];
  desca[ctxt_] = -1;
  BlacsProcessGrid localGrid = BlacsProcessGrid(1, 1);
  if (grid->getCurrentProcess() == 0) {
    desca[ctxt_] = localGrid.getContextHandle();
    desca[m_] = globalRows;
    desca[n_] = globalColumns;
    desca[mb_] = globalRows;
    desca[nb_] = globalColumns;
    desca[rsrc_] = 0;
    desca[csrc_] = 0;
    desca[lld_] = globalRows;
  }
  pdgemr2d_(globalRows, globalColumns, this->getLocalPointer(), 1, 1, this->descriptor,
            localMatrix.getPointer(), 1, 1, desca, this->grid->getContextHandle());

  return localMatrix;
}

double* DataMatrixDistributed::getLocalPointer() { return this->localData.data(); }

const double* DataMatrixDistributed::getLocalPointer() const { return this->localData.data(); }

int* DataMatrixDistributed::getDescriptor() { return this->descriptor; }

const int* DataMatrixDistributed::getDescriptor() const { return this->descriptor; }

int DataMatrixDistributed::getLocalRows() const { return this->localRows; }

int DataMatrixDistributed::getLocalColumns() const { return this->localColumns; }

const double& DataMatrixDistributed::operator()(size_t row, size_t col) const {
  return this->get(row, col);
}

double DataMatrixDistributed::get(size_t row, size_t col) const {
  size_t processRow = globalToProcessIndex(row, grid->getTotalRows(), rowBlockSize, 0);
  size_t processColumn = globalToProcessIndex(row, grid->getTotalColumns(), columnBlockSize, 0);

  double value;
  if (grid->getCurrentColumn() == processColumn && grid->getCurrentRow() == processRow) {
    size_t localRowIndex = globalToLocalIndex(row, grid->getTotalRows(), rowBlockSize);
    size_t localColumnIndex = globalToLocalIndex(col, grid->getTotalColumns(), columnBlockSize);

    value = getLocalPointer()[(localRowIndex * localColumns) + localColumnIndex];

    // broadcast value to other processes so that every process returns the same value
    dgebs2d_(grid->getContextHandle(), "All", "T", 1, 1, &value, 1);
  } else {
    dgebr2d_(grid->getContextHandle(), "All", "T", 1, 1, &value, 1, 0, 0);
  }

  return value;
}

void DataMatrixDistributed::set(size_t row, size_t col, double value) {
  size_t processRow = globalToProcessIndex(row, grid->getTotalRows(), rowBlockSize, 0);
  size_t processColumn = globalToProcessIndex(row, grid->getTotalColumns(), columnBlockSize, 0);

  if (grid->getCurrentColumn() == processColumn && grid->getCurrentRow() == processRow) {
    size_t localRowIndex = globalToLocalIndex(row, grid->getTotalRows(), rowBlockSize);
    size_t localColumnIndex = globalToLocalIndex(col, grid->getTotalColumns(), columnBlockSize);

    getLocalPointer()[(localRowIndex * localColumns) + localColumnIndex] = value;
  }
}

void DataMatrixDistributed::mult(const DataVectorDistributed& x, DataVectorDistributed& y,
                                 double alpha, double beta) const {
  // transpose by default, as pdgemv uses fortran memory order
  // pdgemv: sub(y) := scalar*sub(A)'*sub(x) + beta*sub(y);
  pdgemv_(pblasTranspose, globalRows, globalColumns, 1.0, x.getLocalPointer(), 1, 1,
          x.getDescriptor(), x.getLocalPointer(), 1, 1, x.getDescriptor(), 1, beta,
          y.getLocalPointer(), 1, 1, y.getDescriptor(), 1);
  // TODO(jan) incx, incy 1 oder m_x?
}

void DataMatrixDistributed::mult(bool transposeA, bool transposeB, const DataMatrixDistributed& b,
                                 DataMatrixDistributed& c, double alpha, double beta) const {
  // sub(C) := alpha*op(sub(A))*op(sub(B)) + beta*sub(C)

  // reverse transpose to account for fortran arrays
  const char* transA = (transposeA ? pblasNoTranspose : pblasTranspose);
  const char* transB = (transposeB ? pblasNoTranspose : pblasTranspose);

  pdgemm_(transA, transB, this->getGlobalRows(), b.getGlobalCols(), this->getGlobalCols(), alpha,
          this->getLocalPointer(), 1, 1, this->getDescriptor(), b.getLocalPointer(), 1, 1,
          b.getDescriptor(), beta, c.getLocalPointer(), 1, 1, c.getDescriptor());
}

void DataMatrixDistributed::printMatrix(size_t processRow, size_t processColumn) const {
  std::vector<double> work(globalRows * globalColumns);
  pdlaprnt_(globalRows, globalColumns, getLocalPointer(), 1, 1, getDescriptor(), processRow,
            processColumn, "", 0, work.data());
}

size_t DataMatrixDistributed::globalToLocalIndex(size_t globalIndex, size_t numberOfProcesses,
                                                 size_t blockSize) const {
  // note that the division is rounded
  size_t localBlockIndex = globalIndex / (numberOfProcesses * blockSize);
  size_t elementIndex = globalIndex % blockSize;
  return (localBlockIndex * blockSize) + elementIndex;
}

size_t DataMatrixDistributed::globalToProcessIndex(size_t globalIndex, size_t numberOfProcesses,
                                                   size_t blockSize, size_t processOffset) const {
  size_t processIndex = globalIndex / blockSize;
  return (processOffset + processIndex) % numberOfProcesses;
}

}  // namespace datadriven
}  // namespace sgpp

#endif  // USE_SCALAPACK