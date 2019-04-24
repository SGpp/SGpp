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
#include <sgpp/datadriven/scalapack/DataMatrixDistributed.hpp>

#include <algorithm>
#include <iostream>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>
#include <sgpp/datadriven/scalapack/blacs.hpp>

namespace sgpp {
namespace datadriven {

DataMatrixDistributed::DataMatrixDistributed() {
#ifndef USE_SCALAPACK
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* not USE_SCALAPACK */
}

DataMatrixDistributed::DataMatrixDistributed(std::shared_ptr<BlacsProcessGrid> grid,
                                             size_t globalRows, size_t globalColumns,
                                             size_t rowBlockSize, size_t columnBlockSize,
                                             double value, DataMatrixDistributed::DTYPE dtype)
    : grid(grid),
      globalRows(globalColumns),  // transpose for column-major layout
      globalColumns(globalRows),  // transpose for column-major layour
      rowBlockSize(columnBlockSize),
      columnBlockSize(rowBlockSize),
      localRows(0),
      localColumns(0),
      leadingDimension(1) {
#ifdef USE_SCALAPACK
  descriptor[dtype_] = static_cast<int>(dtype);

  if (isProcessMapped()) {
    if (this->rowBlockSize == 0 || this->columnBlockSize == 0) {
      throw sgpp::base::algorithm_exception("DataMatrixDistributed: block size has to be > 0");
    }

    localRows = numroc_(this->globalRows, this->rowBlockSize, grid->getCurrentRow(), 0,
                        grid->getTotalRows());
    localColumns = numroc_(this->globalColumns, this->columnBlockSize, grid->getCurrentColumn(), 0,
                           grid->getTotalColumns());

    if (localRows == 0 || localColumns == 0) {
      localRows = 0;
      localColumns = 0;
    }

    leadingDimension = std::max<size_t>(1, localRows);

    // initialize local matrix
    localData.assign(localRows * localColumns, value);

    descriptor[ctxt_] = grid->getContextHandle();
    descriptor[m_] = static_cast<int>(this->globalRows);
    descriptor[n_] = static_cast<int>(this->globalColumns);
    descriptor[mb_] = static_cast<int>(this->rowBlockSize);
    descriptor[nb_] = static_cast<int>(this->columnBlockSize);
    descriptor[rsrc_] = 0;
    descriptor[csrc_] = 0;
    descriptor[lld_] = static_cast<int>(leadingDimension);
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

DataMatrixDistributed::DataMatrixDistributed(double* input, std::shared_ptr<BlacsProcessGrid> grid,
                                             size_t globalRows, size_t globalColumns,
                                             size_t rowBlockSize, size_t columnBlockSize,
                                             DataMatrixDistributed::DTYPE dtype)
    : DataMatrixDistributed(grid, globalRows, globalColumns, columnBlockSize, rowBlockSize, 0.0,
                            dtype) {
  if (isProcessMapped()) {
    distribute(input);
  }
}

DataMatrixDistributed DataMatrixDistributed::fromSharedData(const double* input,
                                                            std::shared_ptr<BlacsProcessGrid> grid,
                                                            size_t globalRows, size_t globalColumns,
                                                            size_t rowBlockSize,
                                                            size_t columnBlockSize, DTYPE dtype) {
  DataMatrixDistributed matrix(grid, globalRows, globalColumns, rowBlockSize, columnBlockSize, 0.0,
                               dtype);
  for (size_t row = 0; row < matrix.globalRows; row += matrix.rowBlockSize) {
    for (size_t col = 0; col < matrix.globalColumns; col += matrix.columnBlockSize) {
      // last block might be smaller
      size_t rowsToSend = matrix.rowBlockSize;
      if (row + rowsToSend > matrix.globalRows) {
        rowsToSend = matrix.globalRows - row;
      }

      size_t colsToSend = matrix.columnBlockSize;
      if (col + colsToSend > matrix.globalColumns) {
        colsToSend = matrix.globalColumns - col;
      }

      int receiverRow =
          matrix.globalToProcessIndex(row, grid->getTotalRows(), matrix.rowBlockSize, 0);
      int receiverCol =
          matrix.globalToProcessIndex(col, grid->getTotalColumns(), matrix.columnBlockSize, 0);

      // copy the data to the local matrix of the right process
      if (grid->getCurrentRow() == receiverRow && grid->getCurrentColumn() == receiverCol) {
        size_t localRowIndex =
            matrix.globalToLocalIndex(row, grid->getTotalRows(), matrix.rowBlockSize);
        size_t localColumnIndex =
            matrix.globalToLocalIndex(col, grid->getTotalColumns(), matrix.columnBlockSize);

        for (size_t i = 0; i < colsToSend; i++) {
          const double* submatrix = &input[((col + i) * matrix.globalRows) + row];

          double* localBlock =
              &matrix
                   .getLocalPointer()[((localColumnIndex + i) * matrix.localRows) + localRowIndex];

          std::copy_n(submatrix, rowsToSend, localBlock);
        }
      }
    }
  }
  return matrix;
}

DataMatrix DataMatrixDistributed::toLocalDataMatrix() const {
  DataMatrix localMatrix;
  if (isProcessMapped()) {
    localMatrix = gather();
  }
  return localMatrix;
}

DataMatrix DataMatrixDistributed::toLocalDataMatrixBroadcast() const {
  DataMatrix localMatrix;
  if (isProcessMapped()) {
    localMatrix = broadcast();
  }
  return localMatrix;
}

void DataMatrixDistributed::toLocalDataMatrix(DataMatrix& localMatrix) const {
  if (isProcessMapped()) {
    gather(localMatrix);
  }
}

void DataMatrixDistributed::toLocalDataMatrixBroadcast(DataMatrix& localMatrix) const {
  if (isProcessMapped()) {
    broadcast(localMatrix);
  }
}

double* DataMatrixDistributed::getLocalPointer() { return this->localData.data(); }

const double* DataMatrixDistributed::getLocalPointer() const { return this->localData.data(); }

int* DataMatrixDistributed::getDescriptor() { return this->descriptor; }

const int* DataMatrixDistributed::getDescriptor() const { return this->descriptor; }

size_t DataMatrixDistributed::getLocalRows() const {
  // internal storage is transposed
  return this->localColumns;
}

size_t DataMatrixDistributed::getLocalColumns() const {
  // internal storage is transposed
  return this->localRows;
}

size_t DataMatrixDistributed::getRowBlockSize() const { return this->rowBlockSize; }

size_t DataMatrixDistributed::getColumnBlockSize() const { return this->columnBlockSize; }

size_t DataMatrixDistributed::getGlobalRows() const {
  // internal storage is transposed
  return globalColumns;
}

size_t DataMatrixDistributed::getGlobalCols() const {
  // internal storage is transposed
  return globalRows;
}

std::shared_ptr<BlacsProcessGrid> DataMatrixDistributed::getProcessGrid() const {
  return this->grid;
}

double DataMatrixDistributed::operator()(size_t row, size_t col) const {
  return this->get(row, col);
}

double DataMatrixDistributed::get(size_t row, size_t col) const {
#ifdef USE_SCALAPACK
  // swap row and col to account for transposed storage
  int processRow = globalToProcessIndex(col, grid->getTotalRows(), rowBlockSize, 0);
  int processColumn = globalToProcessIndex(row, grid->getTotalColumns(), columnBlockSize, 0);

  double value = 0.0;
  if (grid->getCurrentColumn() == processColumn && grid->getCurrentRow() == processRow) {
    size_t localRowIndex = globalToLocalIndex(col, grid->getTotalRows(), rowBlockSize);
    size_t localColumnIndex = globalToLocalIndex(row, grid->getTotalColumns(), columnBlockSize);

    value = getLocalPointer()[(localColumnIndex * localRows) + localRowIndex];

    // broadcast value to other processes so that every process returns the same value
    Cdgebs2d(grid->getContextHandle(), "All", "T", 1, 1, &value, 1);
  } else if (isProcessMapped()) {
    Cdgebr2d(grid->getContextHandle(), "All", "T", 1, 1, &value, 1, processRow, processColumn);
  } else {
    std::cout
        << "Warning! Process not in the grid tried to call get, invalid result will be returned!"
        << std::endl;
  }

  return value;
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

void DataMatrixDistributed::set(size_t row, size_t col, double value) {
  // swap row and col to account for transposed storage
  int processRow = globalToProcessIndex(col, grid->getTotalRows(), rowBlockSize, 0);
  int processColumn = globalToProcessIndex(row, grid->getTotalColumns(), columnBlockSize, 0);

  if (grid->getCurrentColumn() == processColumn && grid->getCurrentRow() == processRow) {
    size_t localRowIndex = globalToLocalIndex(col, grid->getTotalRows(), rowBlockSize);
    size_t localColumnIndex = globalToLocalIndex(row, grid->getTotalColumns(), columnBlockSize);

    getLocalPointer()[(localColumnIndex * localRows) + localRowIndex] = value;
  }
}

void DataMatrixDistributed::setAll(double value) {
  if (isProcessMapped()) {
    std::fill(localData.begin(), localData.end(), value);
  }
}

void DataMatrixDistributed::copyFrom(const DataMatrixDistributed& other) {
  // use std::vector copy constructor
  this->localData = other.localData;

  // update the dimensions and the descriptor
  this->localRows = other.localRows;
  this->localColumns = other.localColumns;

  this->globalRows = other.globalRows;
  this->globalColumns = other.globalColumns;

  this->leadingDimension = other.leadingDimension;

  descriptor[m_] = static_cast<int>(this->globalRows);
  descriptor[n_] = static_cast<int>(this->globalColumns);
  descriptor[lld_] = static_cast<int>(this->leadingDimension);
}

DataMatrixDistributed DataMatrixDistributed::transpose() {
  // note that globalRows/columns and row/columnBlockSize are switched
  DataMatrixDistributed transposed =
      DataMatrixDistributed(grid, globalRows, globalColumns, rowBlockSize, columnBlockSize);
  DataMatrixDistributed::transpose(*this, transposed);
  return transposed;
}

void DataMatrixDistributed::transpose(const DataMatrixDistributed& a, DataMatrixDistributed& c,
                                      double alpha, double beta) {
#ifdef USE_SCALAPACK
  if (a.isProcessMapped() || c.isProcessMapped()) {
    pdtran_(c.getGlobalCols(), c.getGlobalRows(), alpha, a.getLocalPointer(), 1, 1,
            a.getDescriptor(), beta, c.getLocalPointer(), 1, 1, c.getDescriptor());
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

void DataMatrixDistributed::add(const DataMatrixDistributed& a) {
  DataMatrixDistributed::add(*this, a);
}

void DataMatrixDistributed::add(DataMatrixDistributed& c, const DataMatrixDistributed& a,
                                bool transposeA, double beta, double alpha) {
#ifdef USE_SCALAPACK
  if (a.isProcessMapped() || c.isProcessMapped()) {
    const char* transA = (transposeA ? pblasTranspose : pblasNoTranspose);

    pdgeadd_(transA, a.getGlobalCols(), a.getGlobalRows(), alpha, a.getLocalPointer(), 1, 1,
             a.getDescriptor(), beta, c.getLocalPointer(), 1, 1, c.getDescriptor());
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

void DataMatrixDistributed::sub(const DataMatrixDistributed& a) {
  DataMatrixDistributed::sub(*this, a);
}

void DataMatrixDistributed::sub(DataMatrixDistributed& c, const DataMatrixDistributed& a,
                                bool transa, double beta, double alpha) {
  DataMatrixDistributed::add(c, a, transa, beta, -alpha);
}

void DataMatrixDistributed::mult(const DataVectorDistributed& x, DataVectorDistributed& y,
                                 bool transpose, double alpha, double beta) const {
  DataMatrixDistributed::mult(*this, x, y, transpose, alpha, beta);
}

void DataMatrixDistributed::mult(const DataMatrixDistributed& a, const DataVectorDistributed& x,
                                 DataVectorDistributed& y, bool transpose, double alpha,
                                 double beta) {
#ifdef USE_SCALAPACK
  if (a.isProcessMapped() || x.isProcessMapped() || y.isProcessMapped()) {
    // transpose matrix by default, as internal storage is transposed
    const char* trans = (transpose ? pblasNoTranspose : pblasTranspose);

    // pdgemv: sub(y) := scalar*sub(A)'*sub(x) + beta*sub(y);
    pdgemv_(trans, a.getGlobalCols(), a.getGlobalRows(), 1.0, a.getLocalPointer(), 1, 1,
            a.getDescriptor(), x.getLocalPointer(), 1, 1, x.getDescriptor(), 1, beta,
            y.getLocalPointer(), 1, 1, y.getDescriptor(), 1);
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

void DataMatrixDistributed::mult(const DataMatrixDistributed& b, DataMatrixDistributed& c,
                                 bool transposeA, bool transposeB, double alpha,
                                 double beta) const {
  DataMatrixDistributed::mult(*this, b, c, transposeA, transposeB, alpha, beta);
}

void DataMatrixDistributed::mult(const DataMatrixDistributed& a, const DataMatrixDistributed& b,
                                 DataMatrixDistributed& c, bool transposeA, bool transposeB,
                                 double alpha, double beta) {
#ifdef USE_SCALAPACK
  if (a.isProcessMapped() || b.isProcessMapped() || c.isProcessMapped()) {
    const char* transA = (transposeA ? pblasTranspose : pblasNoTranspose);
    const char* transB = (transposeB ? pblasTranspose : pblasNoTranspose);

    pdgemm_(transB, transA, b.getGlobalCols(), a.getGlobalRows(), b.getGlobalRows(), alpha,
            b.getLocalPointer(), 1, 1, b.getDescriptor(), a.getLocalPointer(), 1, 1,
            a.getDescriptor(), beta, c.getLocalPointer(), 1, 1, c.getDescriptor());
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

void DataMatrixDistributed::solveCholesky(const DataMatrixDistributed& l, DataVectorDistributed& b,
                                          DataMatrixDistributed::TRIANGULAR uplo) {
#ifdef USE_SCALAPACK
  if (l.isProcessMapped() || b.isProcessMapped()) {
    // implementation detail: values of upper and lower are switched, as internally ScaLAPACK
    // transposes all matrices. By switching upper and lower, everything works as expected.
    const char* tri = (uplo == TRIANGULAR::LOWER ? upperTriangular : lowerTriangular);

    int info = 0;
    pdpotrs_(tri, l.getGlobalRows(), 1, l.getLocalPointer(), 1, 1, l.getDescriptor(),
             b.getLocalPointer(), 1, 1, b.getDescriptor(), info);

    if (info < 0) {
      throw sgpp::base::algorithm_exception("DataMatrixDistributed::solveCholesky() failed");
    }
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

void DataMatrixDistributed::resize(size_t rows, size_t cols) {
#ifdef USE_SCALAPACK
  if (getGlobalRows() == rows && getGlobalCols() == cols) {
    return;
  }
  // transpose
  this->globalRows = cols;
  this->globalColumns = rows;

  if (isProcessMapped()) {
    localRows = numroc_(this->globalRows, this->rowBlockSize, grid->getCurrentRow(), 0,
                        grid->getTotalRows());
    localColumns = numroc_(this->globalColumns, this->columnBlockSize, grid->getCurrentColumn(), 0,
                           grid->getTotalColumns());

    if (localRows == 0 || localColumns == 0) {
      localRows = 0;
      localColumns = 0;
    }

    leadingDimension = std::max<size_t>(1, localRows);

    // resize local matrix
    localData.resize(localRows * localColumns);

    descriptor[m_] = static_cast<int>(this->globalRows);
    descriptor[n_] = static_cast<int>(this->globalColumns);
    descriptor[lld_] = static_cast<int>(leadingDimension);
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

DataVectorDistributed DataMatrixDistributed::toVector() const {
  if (globalRows != 1 && globalColumns != 1) {
    throw base::algorithm_exception(
        "matrix with more than one rows and columns can not be converted to a vector");
  }
  DataMatrix localMatrix;
  if (grid->getCurrentRow() == 0 && grid->getCurrentColumn() == 0) {
    localMatrix = this->toLocalDataMatrix();
  }
  size_t globalSize = (globalRows == 1 ? globalColumns : globalRows);
  return DataVectorDistributed(localMatrix.data(), grid, globalSize, rowBlockSize);
}

void DataMatrixDistributed::printMatrix() const {
  DataMatrix localMatrix = toLocalDataMatrix();
  if (grid->getCurrentProcess() == 0) {
    std::cout << localMatrix.toString() << std::endl;
  }
}

bool DataMatrixDistributed::isProcessMapped() const { return grid->isProcessInGrid(); }

size_t DataMatrixDistributed::globalToLocalRowIndex(size_t globalRowIndex) const {
  return globalToLocalIndex(globalRowIndex, grid->getTotalColumns(), columnBlockSize);
}

size_t DataMatrixDistributed::globalToLocalColumnIndex(size_t globalColumnIndex) const {
  return globalToLocalIndex(globalColumnIndex, grid->getTotalRows(), rowBlockSize);
}

size_t DataMatrixDistributed::localToGlobalRowIndex(size_t localRowIndex) const {
  return localToGlobalIndex(localRowIndex, grid->getCurrentColumn(), grid->getTotalColumns(),
                            columnBlockSize);
}

size_t DataMatrixDistributed::localToGlobalColumnIndex(size_t localColumnIndex) const {
  return localToGlobalIndex(localColumnIndex, grid->getCurrentRow(), grid->getTotalRows(),
                            rowBlockSize);
}

size_t DataMatrixDistributed::globalToRowProcessIndex(size_t globalRowIndex) const {
  return globalToProcessIndex(globalRowIndex, grid->getTotalColumns(), columnBlockSize);
}

size_t DataMatrixDistributed::globalToColumnProcessIndex(size_t globalColumnIndex) const {
  return globalToProcessIndex(globalColumnIndex, grid->getTotalRows(), rowBlockSize);
}

void DataMatrixDistributed::distribute(const double* matrix, int masterRow, int masterCol) {
#ifdef USE_SCALAPACK
  for (size_t row = 0; row < globalRows; row += rowBlockSize) {
    for (size_t col = 0; col < globalColumns; col += columnBlockSize) {
      // last block might be smaller
      size_t rowsToSend = rowBlockSize;
      if (row + rowsToSend > globalRows) {
        rowsToSend = globalRows - row;
      }

      size_t colsToSend = columnBlockSize;
      if (col + colsToSend > globalColumns) {
        colsToSend = globalColumns - col;
      }

      int receiverRow = globalToProcessIndex(row, grid->getTotalRows(), rowBlockSize, 0);
      int receiverCol = globalToProcessIndex(col, grid->getTotalColumns(), columnBlockSize, 0);

      // send before receive to prevent deadlocks
      if (grid->getCurrentRow() == masterRow && grid->getCurrentColumn() == masterCol) {
        const double* submatrix = &matrix[(col * globalRows) + row];

        Cdgesd2d(grid->getContextHandle(), rowsToSend, colsToSend, submatrix, globalRows,
                 receiverRow, receiverCol);
      }

      // receive (sends and receives to the same process are strictly ordered, so indexing works)
      if (grid->getCurrentRow() == receiverRow && grid->getCurrentColumn() == receiverCol) {
        size_t localRowIndex = globalToLocalIndex(row, grid->getTotalRows(), rowBlockSize);
        size_t localColumnIndex = globalToLocalIndex(col, grid->getTotalColumns(), columnBlockSize);

        double* localBlock = &getLocalPointer()[(localColumnIndex * localRows) + localRowIndex];

        // switch between rows and columns to account for column-major memory access
        Cdgerv2d(grid->getContextHandle(), rowsToSend, colsToSend, localBlock, leadingDimension,
                 masterRow, masterCol);
      }
    }
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

void DataMatrixDistributed::gather(DataMatrix& localMatrix, int masterRow, int masterCol) const {
#ifdef USE_SCALAPACK
  for (size_t row = 0; row < globalRows; row += rowBlockSize) {
    for (size_t col = 0; col < globalColumns; col += columnBlockSize) {
      // last block might be smaller
      size_t rowsToSend = rowBlockSize;
      if (row + rowsToSend > globalRows) {
        rowsToSend = globalRows - row;
      }

      size_t colsToSend = columnBlockSize;
      if (col + colsToSend > globalColumns) {
        colsToSend = globalColumns - col;
      }

      int senderRow = globalToProcessIndex(row, grid->getTotalRows(), rowBlockSize, 0);
      int senderCol = globalToProcessIndex(col, grid->getTotalColumns(), columnBlockSize, 0);

      // send before receive to prevent deadlocks
      if (grid->getCurrentRow() == senderRow && grid->getCurrentColumn() == senderCol) {
        size_t localRowIndex = globalToLocalIndex(row, grid->getTotalRows(), rowBlockSize);
        size_t localColumnIndex = globalToLocalIndex(col, grid->getTotalColumns(), columnBlockSize);

        const double* localBlock =
            &getLocalPointer()[(localColumnIndex * localRows) + localRowIndex];

        // switch between rows and columns to account for column-major memory access
        Cdgesd2d(grid->getContextHandle(), rowsToSend, colsToSend, localBlock, leadingDimension,
                 masterRow, masterCol);
      }

      // receive (sends and receives to the same process are strictly ordered, so indexing works)
      if (grid->getCurrentRow() == masterRow && grid->getCurrentColumn() == masterCol) {
        double* submatrix = &localMatrix.getPointer()[(col * globalRows) + row];

        Cdgerv2d(grid->getContextHandle(), rowsToSend, colsToSend, submatrix, globalRows, senderRow,
                 senderCol);
      }
    }
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

DataMatrix DataMatrixDistributed::gather(int masterRow, int masterCol) const {
  DataMatrix localMatrix;
  if (grid->getCurrentRow() == masterRow && grid->getCurrentColumn() == masterCol) {
    // only init localMatrix on master process
    localMatrix = DataMatrix(globalColumns, globalRows);  // transpose
  }
  gather(localMatrix, masterRow, masterCol);
  return localMatrix;
}

void DataMatrixDistributed::broadcast(DataMatrix& localMatrix) const {
#ifdef USE_SCALAPACK
  for (size_t row = 0; row < globalRows; row += rowBlockSize) {
    for (size_t col = 0; col < globalColumns; col += columnBlockSize) {
      // last block might be smaller
      size_t rowsToSend = rowBlockSize;
      if (row + rowsToSend > globalRows) {
        rowsToSend = globalRows - row;
      }

      size_t colsToSend = columnBlockSize;
      if (col + colsToSend > globalColumns) {
        colsToSend = globalColumns - col;
      }

      int senderRow = globalToProcessIndex(row, grid->getTotalRows(), rowBlockSize, 0);
      int senderCol = globalToProcessIndex(col, grid->getTotalColumns(), columnBlockSize, 0);

      // send before receive to prevent deadlocks
      if (grid->getCurrentRow() == senderRow && grid->getCurrentColumn() == senderCol) {
        size_t localRowIndex = globalToLocalIndex(row, grid->getTotalRows(), rowBlockSize);
        size_t localColumnIndex = globalToLocalIndex(col, grid->getTotalColumns(), columnBlockSize);

        const double* localBlock =
            &getLocalPointer()[(localColumnIndex * localRows) + localRowIndex];

        // switch between rows and columns to account for column-major memory access
        Cdgebs2d(grid->getContextHandle(), "All", "T", rowsToSend, colsToSend, localBlock,
                 leadingDimension);

        // for the sender, directly copy the block into the local matrix
        for (size_t i = 0; i < colsToSend; i++) {
          double* submatrix = &localMatrix.getPointer()[((col + i) * globalRows) + row];
          std::copy_n(localBlock, rowsToSend * colsToSend, submatrix);

          localBlock += localRows;
        }
      } else {
        // receive on all other processes (sends and receives to the same process are strictly
        // ordered, so indexing works)
        double* submatrix = &localMatrix.getPointer()[(col * globalRows) + row];

        Cdgebr2d(grid->getContextHandle(), "All", "T", rowsToSend, colsToSend, submatrix,
                 globalRows, senderRow, senderCol);
      }
    }
  }
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

DataMatrix DataMatrixDistributed::broadcast() const {
  DataMatrix localMatrix = DataMatrix(globalColumns, globalRows);  // transpose
  broadcast(localMatrix);
  return localMatrix;
}

size_t DataMatrixDistributed::globalToLocalIndex(size_t globalIndex, size_t numberOfProcesses,
                                                 size_t blockSize) const {
  // note that the division is rounded down
  size_t localBlockIndex = globalIndex / (numberOfProcesses * blockSize);
  size_t elementIndex = globalIndex % blockSize;
  return (localBlockIndex * blockSize) + elementIndex;
}

size_t DataMatrixDistributed::localToGlobalIndex(size_t localIndex, size_t process,
                                                 size_t numberOfProcesses, size_t blockSize) const {
  // note that the division is rounded down
  size_t localBlockIndex = localIndex / blockSize;
  size_t elementIndex = localIndex % blockSize;
  return (localBlockIndex * numberOfProcesses * blockSize) + process * blockSize + elementIndex;
}

int DataMatrixDistributed::globalToProcessIndex(size_t globalIndex, size_t numberOfProcesses,
                                                size_t blockSize, int processOffset) const {
  size_t processIndex = globalIndex / blockSize;
  size_t process = (processOffset + processIndex) % numberOfProcesses;
  return static_cast<int>(process);
}

}  // namespace datadriven
}  // namespace sgpp