/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataVectorDistributed.cpp
 *
 * Created on: Feb 19, 2019
 *     Author: Jan Schopohl
 */
#ifdef USE_SCALAPACK

#include <sgpp/datadriven/scalapack/DataVectorDistributed.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataVector;

DataVectorDistributed::DataVectorDistributed(std::shared_ptr<BlacsProcessGrid> grid,
                                             size_t globalSize, size_t blockSize, double value)
    : data(grid, 1, globalSize, blockSize, blockSize, value, DataMatrixDistributed::DTYPE::DENSE),
      grid(grid) {
  // internal storage in DataMatrixDistributed in columns, as fortran uses column-major order
}

DataVectorDistributed::DataVectorDistributed(double* input, std::shared_ptr<BlacsProcessGrid> grid,
                                             size_t globalSize, size_t blockSize)
    : data(input, grid, 1, globalSize, blockSize, blockSize, DataMatrixDistributed::DTYPE::DENSE),
      grid(grid) {
  // internal storage in DataMatrixDistributed in columns, as fortran uses column-major order
}

double DataVectorDistributed::operator()(size_t row) const { return data(0, row); }

double DataVectorDistributed::get(size_t row) const { return data.get(0, row); }

void DataVectorDistributed::set(size_t row, double value) { data.set(0, row, value); }

void DataVectorDistributed::setAll(double value) { data.setAll(value); }

void DataVectorDistributed::add(const DataVectorDistributed& x) {
  DataVectorDistributed::add(*this, x);
}

void DataVectorDistributed::add(DataVectorDistributed& y, const DataVectorDistributed& x,
                                double a) {
  if (y.isProcessMapped() || x.isProcessMapped()) {
    pdaxpy_(y.getGlobalRows(), a, x.getLocalPointer(), 1, 1, x.getDescriptor(), 1,
            y.getLocalPointer(), 1, 1, y.getDescriptor(), 1);
  }
}

double DataVectorDistributed::dot(const DataVectorDistributed& y) const {
  return DataVectorDistributed::dot(*this, y);
}

double DataVectorDistributed::dot(const DataVectorDistributed& x, const DataVectorDistributed& y) {
  double dot = 0.0;
  if (x.isProcessMapped() || y.isProcessMapped()) {
    pddot_(x.getGlobalRows(), dot, x.getLocalPointer(), 1, 1, x.getDescriptor(), 1,
           y.getLocalPointer(), 1, 1, y.getDescriptor(), 1);
  }
  return dot;
}

void DataVectorDistributed::scale(double a) {
  if (isProcessMapped()) {
    pdscal_(getGlobalRows(), a, getLocalPointer(), 1, 1, getDescriptor(), 1);
  }
}

void DataVectorDistributed::append(size_t rows) { data.appendRows(rows); }

void DataVectorDistributed::resize(size_t rows) { data.resize(1, rows); }

double* DataVectorDistributed::getLocalPointer() { return data.getLocalPointer(); }

const double* DataVectorDistributed::getLocalPointer() const { return data.getLocalPointer(); }

std::shared_ptr<BlacsProcessGrid> DataVectorDistributed::getProcessGrid() const {
  return data.getProcessGrid();
}

DataVector DataVectorDistributed::toLocalDataVector() const {
  DataMatrix localMatrix = data.toLocalDataMatrix();
  if (grid->getCurrentRow() == 0 && grid->getCurrentColumn() == 0) {
    DataVector localVector(localMatrix.data(), localMatrix.getNcols());
    return localVector;
  }
  return DataVector();
}

DataVector DataVectorDistributed::toLocalDataVectorBroadcast() const {
  DataMatrix localMatrix = data.toLocalDataMatrixBroadcast();
  if (grid->isProcessInGrid()) {
    DataVector localVector(localMatrix.data(), localMatrix.getNcols());
    return localVector;
  }
  return DataVector();
}

void DataVectorDistributed::toLocalDataVector(DataVector& localVector) const {
  DataMatrix localMatrix = data.toLocalDataMatrix();
  if (grid->getCurrentRow() == 0 && grid->getCurrentColumn() == 0) {
    localVector.assign(localMatrix.begin(), localMatrix.end());
  }
}

void DataVectorDistributed::toLocalDataVectorBroadcast(DataVector& localVector) const {
  DataMatrix localMatrix = data.toLocalDataMatrixBroadcast();
  if (grid->isProcessInGrid()) {
    localVector.assign(localMatrix.begin(), localMatrix.end());
  }
}

int* DataVectorDistributed::getDescriptor() { return data.getDescriptor(); }

const int* DataVectorDistributed::getDescriptor() const { return data.getDescriptor(); }

size_t DataVectorDistributed::getGlobalRows() const { return data.getGlobalCols(); }

size_t DataVectorDistributed::getLocalRows() const { return data.getLocalColumns(); }

size_t DataVectorDistributed::getBlockSize() const { return data.getColumnBlockSize(); }

void DataVectorDistributed::printVector() const {
  DataVector localVector = toLocalDataVector();
  if (grid->getCurrentRow() == 0 && grid->getCurrentColumn() == 0) {
    std::cout << localVector.toString() << std::endl;
  }
}

bool DataVectorDistributed::isProcessMapped() const { return data.isProcessMapped(); }

DataMatrixDistributed& DataVectorDistributed::getMatrix() { return data; }

const DataMatrixDistributed& DataVectorDistributed::getMatrix() const { return data; }

void DataVectorDistributed::distribute(double* input, int masterRow, int masterCol) {
  data.distribute(input, masterRow, masterCol);
}

size_t DataVectorDistributed::globalToLocalRowIndex(size_t globalRowIndex) const {
  return data.globalToLocalColumnIndex(globalRowIndex);
}

size_t DataVectorDistributed::localToGlobalRowIndex(size_t localRowIndex) const {
  return data.localToGlobalColumnIndex(localRowIndex);
}

size_t DataVectorDistributed::globalToRowProcessIndex(size_t globalRowIndex) const {
  return data.globalToColumnProcessIndex(globalRowIndex);
}

}  // namespace datadriven
}  // namespace sgpp

#endif