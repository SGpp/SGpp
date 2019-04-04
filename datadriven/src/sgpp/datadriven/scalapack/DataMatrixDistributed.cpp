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

#include <algorithm>
#include <iostream>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>

#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>
#include <sgpp/datadriven/scalapack/blacs.hpp>

namespace sgpp {
namespace datadriven {

DataMatrixDistributed::DataMatrixDistributed() {}

DataMatrixDistributed::DataMatrixDistributed(std::shared_ptr<BlacsProcessGrid> grid, int globalRows,
                                             int globalColumns, int rowBlockSize,
                                             int columnBlockSize, double value,
                                             DataMatrixDistributed::DTYPE dtype)
    : grid(grid),
      globalRows(globalColumns),  // transpose for column-major layout
      globalColumns(globalRows),  // transpose for column-major layour
      rowBlockSize(columnBlockSize),
      columnBlockSize(rowBlockSize),
      localRows(0),
      localColumns(0),
      leadingDimension(1) {
  descriptor[dtype_] = static_cast<int>(dtype);

  if (isProcessMapped()) {
    localRows = numroc_(this->globalRows, this->rowBlockSize, grid->getCurrentRow(), 0,
                        grid->getTotalRows());
    localColumns = numroc_(this->globalColumns, this->columnBlockSize, grid->getCurrentColumn(), 0,
                           grid->getTotalColumns());

    if (localRows == 0 || localColumns == 0) {
      localRows = 0;
      localColumns = 0;
    }

    leadingDimension = std::max(1, localRows);

    // initialize local matrix
    localData.assign(localRows * localColumns, value);

    descriptor[ctxt_] = grid->getContextHandle();
    descriptor[m_] = this->globalRows;
    descriptor[n_] = this->globalColumns;
    descriptor[mb_] = this->rowBlockSize;
    descriptor[nb_] = this->columnBlockSize;
    descriptor[rsrc_] = 0;
    descriptor[csrc_] = 0;
    descriptor[lld_] = leadingDimension;
  }
}

DataMatrixDistributed::DataMatrixDistributed(double* input, std::shared_ptr<BlacsProcessGrid> grid,
                                             int globalRows, int globalColumns, int rowBlockSize,
                                             int columnBlockSize,
                                             DataMatrixDistributed::DTYPE dtype)
    : DataMatrixDistributed(grid, globalRows, globalColumns, columnBlockSize, rowBlockSize, 0.0,
                            dtype) {
  if (isProcessMapped()) {
    distribute(input);
  }
}

DataMatrixDistributed DataMatrixDistributed::fromSharedData(const double* input,
                                                            std::shared_ptr<BlacsProcessGrid> grid,
                                                            int globalRows, int globalColumns,
                                                            int rowBlockSize, int columnBlockSize,
                                                            DTYPE dtype) {
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

        const double* submatrix = &input[(col * matrix.globalRows) + row];
        double* localBlock =
            &matrix.getLocalPointer()[(localColumnIndex * matrix.localRows) + localRowIndex];

        std::copy_n(submatrix, rowsToSend * colsToSend, localBlock);
      }
    }
  }
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

  descriptor[m_] = this->globalRows;
  descriptor[n_] = this->globalColumns;
  descriptor[lld_] = this->leadingDimension;
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
  if (a.isProcessMapped() || c.isProcessMapped()) {
    pdtran_(c.getGlobalCols(), c.getGlobalRows(), alpha, a.getLocalPointer(), 1, 1,
            a.getDescriptor(), beta, c.getLocalPointer(), 1, 1, c.getDescriptor());
  }
}

void DataMatrixDistributed::add(const DataMatrixDistributed& a) {
  DataMatrixDistributed::add(*this, a);
}

void DataMatrixDistributed::add(DataMatrixDistributed& c, const DataMatrixDistributed& a,
                                bool transposeA, double beta, double alpha) {
  if (a.isProcessMapped() || c.isProcessMapped()) {
    const char* transA = (transposeA ? pblasTranspose : pblasNoTranspose);

    pdgeadd_(transA, a.getGlobalCols(), a.getGlobalRows(), alpha, a.getLocalPointer(), 1, 1,
             a.getDescriptor(), beta, c.getLocalPointer(), 1, 1, c.getDescriptor());
  }
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
  if (a.isProcessMapped() || x.isProcessMapped() || y.isProcessMapped()) {
    // transpose matrix by default, as internal storage is transposed
    const char* trans = (transpose ? pblasNoTranspose : pblasTranspose);

    // pdgemv: sub(y) := scalar*sub(A)'*sub(x) + beta*sub(y);
    pdgemv_(trans, a.getGlobalCols(), a.getGlobalRows(), 1.0, a.getLocalPointer(), 1, 1,
            a.getDescriptor(), x.getLocalPointer(), 1, 1, x.getDescriptor(), 1, beta,
            y.getLocalPointer(), 1, 1, y.getDescriptor(), 1);
  }
}

void DataMatrixDistributed::mult(const DataMatrixDistributed& b, DataMatrixDistributed& c,
                                 bool transposeA, bool transposeB, double alpha,
                                 double beta) const {
  DataMatrixDistributed::mult(*this, b, c, transposeA, transposeB, alpha, beta);
}

void DataMatrixDistributed::mult(const DataMatrixDistributed& a, const DataMatrixDistributed& b,
                                 DataMatrixDistributed& c, bool transposeA, bool transposeB,
                                 double alpha, double beta) {
  // sub(C) := alpha*op(sub(A))*op(sub(B)) + beta*sub(C)

  if (a.isProcessMapped() || b.isProcessMapped() || c.isProcessMapped()) {
    const char* transA = (transposeA ? pblasTranspose : pblasNoTranspose);
    const char* transB = (transposeB ? pblasTranspose : pblasNoTranspose);

    // TODO(jan) assert same block size for all matrices
    // TODO(jan) assert that (b.cols == c.cols) and (a.rows == c.rows) and (b.cols = a.rows)
    pdgemm_(transB, transA, b.getGlobalCols(), a.getGlobalRows(), b.getGlobalRows(), alpha,
            b.getLocalPointer(), 1, 1, b.getDescriptor(), a.getLocalPointer(), 1, 1,
            a.getDescriptor(), beta, c.getLocalPointer(), 1, 1, c.getDescriptor());
  }
}

void DataMatrixDistributed::solveCholesky(const DataMatrixDistributed& l, DataVectorDistributed& b,
                                          DataMatrixDistributed::TRIANGULAR uplo) {
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
}

/*void DataMatrixDistributed::appendRows(size_t rows) {
  if (rows == 0) {
    return;
  }
  if (this->localColumns > 1) {
    throw sgpp::base::not_implemented_exception(
        "DataMatrixDistributed::resizeRows currently only works if there is just one localColumn");
  }

  size_t newRows = globalRows + rows;

  for (size_t row = globalRows; row < newRows; row += rowBlockSize) {
    for (size_t col = 0; col < globalColumns; col += columnBlockSize) {
      // first and last block might be smaller

      size_t rowsToInsert = rowBlockSize - (row % rowBlockSize);
      if (row + rowsToInsert > globalRows) {
        rowsToInsert = globalRows - row;
      }

      size_t colsToInsert = columnBlockSize;
      if (col + colsToInsert > globalColumns) {
        colsToInsert = globalColumns - col;
      }

      int processRow = globalToProcessIndex(row, grid->getTotalRows(), rowBlockSize, 0);
      int processCol = globalToProcessIndex(col, grid->getTotalColumns(), columnBlockSize, 0);

      if (grid->getCurrentRow() == processRow && grid->getCurrentColumn() == processCol) {
        localData.insert(localData.end(), colsToInsert * rowsToInsert, 0.0);
      }

      if (row == globalRows) {
        // align to block size if needed
        row -= globalRows % rowBlockSize;
      }
    }
  }

  globalRows = newRows;
}

void DataMatrixDistributed::appendColumns(size_t cols) {
  if (cols == 0) {
    return;
  }

  size_t newColumns = globalColumns + cols;

  // add new columns
  for (size_t row = 0; row < globalRows; row += rowBlockSize) {
    for (size_t col = globalColumns; col < newColumns; col += columnBlockSize) {
      // first and last block might be smaller

      size_t rowsToInsert = rowBlockSize;
      if (row + rowsToInsert > globalRows) {
        rowsToInsert = globalRows - row;
      }

      size_t colsToInsert = columnBlockSize - (col % columnBlockSize);
      if (col + colsToInsert > newColumns) {
        colsToInsert = globalColumns - col;
      }

      int processRow = globalToProcessIndex(row, grid->getTotalRows(), rowBlockSize, 0);
      int processCol = globalToProcessIndex(col, grid->getTotalColumns(), columnBlockSize, 0);

      if (grid->getCurrentRow() == processRow && grid->getCurrentColumn() == processCol) {
        localData.insert(localData.end(), colsToInsert * rowsToInsert, 0.0);
      }

      if (col == globalColumns) {
        // align to block size if needed
        col -= globalColumns % columnBlockSize;
      }
    }
  }

  globalColumns = newColumns;
}*/

void DataMatrixDistributed::resize(size_t rows, size_t cols) {
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

    leadingDimension = std::max(1, localRows);

    // resize local matrix
    localData.resize(localRows * localColumns);

    descriptor[m_] = this->globalRows;
    descriptor[n_] = this->globalColumns;
    descriptor[lld_] = leadingDimension;
  }
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
  return globalToLocalIndex(globalRowIndex, grid->getTotalRows(), rowBlockSize);
}

size_t DataMatrixDistributed::globalToLocalColumnIndex(size_t globalColumnIndex) const {
  return globalToLocalIndex(globalColumnIndex, grid->getTotalColumns(), columnBlockSize);
}

size_t DataMatrixDistributed::localToGlobalRowIndex(size_t localRowIndex) const {
  return localToGlobalIndex(localRowIndex, grid->getCurrentRow(), grid->getTotalRows(),
                            rowBlockSize);
}

size_t DataMatrixDistributed::localToGlobalColumnIndex(size_t localColumnIndex) const {
  return localToGlobalIndex(localColumnIndex, grid->getCurrentColumn(), grid->getTotalColumns(),
                            columnBlockSize);
}

size_t DataMatrixDistributed::globalToRowProcessIndex(size_t globalRowIndex) const {
  return globalToProcessIndex(globalRowIndex, grid->getTotalRows(), rowBlockSize);
}

size_t DataMatrixDistributed::globalToColumnProcessIndex(size_t globalColumnIndex) const {
  return globalToProcessIndex(globalColumnIndex, grid->getTotalColumns(), columnBlockSize);
}

void DataMatrixDistributed::distribute(const double* matrix, int masterRow, int masterCol) {
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
}

void DataMatrixDistributed::gather(DataMatrix& localMatrix, int masterRow, int masterCol) const {
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

        double* submatrix = &localMatrix.getPointer()[(col * globalRows) + row];

        // for the sender, directly copy the block into the local matrix
        std::copy_n(localBlock, rowsToSend * colsToSend, submatrix);
      } else {
        // receive on all other processes (sends and receives to the same process are strictly
        // ordered, so indexing works)
        double* submatrix = &localMatrix.getPointer()[(col * globalRows) + row];

        Cdgebr2d(grid->getContextHandle(), "All", "T", rowsToSend, colsToSend, submatrix,
                 globalRows, senderRow, senderCol);
      }
    }
  }
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

#endif  // USE_SCALAPACK