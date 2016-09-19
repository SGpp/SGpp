// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>
// #include <chrono>
// #include <thread>
#include <iostream>
#include <vector>

#include "OperationMultiEvalHPX.hpp"
#include "sgpp/base/exception/not_implemented_exception.hpp"
#include "sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp"

namespace sgpp {
namespace datadriven {

OperationMultiEvalHPX::OperationMultiEvalHPX(base::Grid& grid,
        base::DataMatrix& dataset, OperationMultipleEvalType nodeImplType,
        OperationMultipleEvalSubType nodeImplSubType, bool verbose) :
        OperationMultipleEval(grid, dataset), nodeImplType(nodeImplType), nodeImplSubType(
                nodeImplSubType), dim(grid.getDimension()), verbose(verbose), duration(
                -1.0) {
    // create the kernel specific data structures for the current grid
    this->prepare();
}

OperationMultiEvalHPX::~OperationMultiEvalHPX() {
}

void OperationMultiEvalHPX::mult(sgpp::base::DataVector& alpha,
        sgpp::base::DataVector& result) {
    throw base::not_implemented_exception();
}

void OperationMultiEvalHPX::multTranspose(sgpp::base::DataVector& source,
        sgpp::base::DataVector& result) {
    throw base::not_implemented_exception();
}

double OperationMultiEvalHPX::getDuration() {
    return this->duration;
}

void OperationMultiEvalHPX::prepare() {
}
}  // namespace datadriven
}  // namespace sgpp
