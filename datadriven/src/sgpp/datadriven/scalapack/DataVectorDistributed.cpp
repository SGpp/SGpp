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

namespace sgpp {
namespace datadriven {

DataVectorDistributed::DataVectorDistributed(std::shared_ptr<BlacsProcessGrid> grid,
                                             size_t globalSize, size_t blockSize, double value)
    : data(grid, globalSize, 1, DataMatrixDistributed::DTYPE::DENSE, blockSize, blockSize, value) {}

double* DataVectorDistributed::getLocalPointer() { return data.getLocalPointer(); }

const double* DataVectorDistributed::getLocalPointer() const { return data.getLocalPointer(); }

int* DataVectorDistributed::getDescriptor() { return data.getDescriptor(); }

const int* DataVectorDistributed::getDescriptor() const { return data.getDescriptor(); }

}  // namespace datadriven
}  // namespace sgpp

#endif