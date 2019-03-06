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

double* DataVectorDistributed::getLocalPointer() { return data.getLocalPointer(); }

const double* DataVectorDistributed::getLocalPointer() const { return data.getLocalPointer(); }

DataVector DataVectorDistributed::toLocalDataVector() const {
  DataMatrix localMatrix = data.toLocalDataMatrix();
  if (grid->getCurrentProcess() == 0) {
    DataVector localVector(localMatrix.data(), localMatrix.getNrows());
    return localVector;
  }
  return DataVector();
}

int* DataVectorDistributed::getDescriptor() { return data.getDescriptor(); }

const int* DataVectorDistributed::getDescriptor() const { return data.getDescriptor(); }

size_t DataVectorDistributed::getGlobalRows() const { return data.getGlobalCols(); }

size_t DataVectorDistributed::getLocalRows() const { return data.getLocalColumns(); }

void DataVectorDistributed::printVector() const {
  DataVector localVector = toLocalDataVector();
  if (grid->getCurrentProcess()) {
    std::cout << localVector.toString() << std::endl;
  }
}

bool DataVectorDistributed::isProcessMapped() const { return data.isProcessMapped(); }

}  // namespace datadriven
}  // namespace sgpp

#endif