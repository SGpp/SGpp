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
#include "sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/OperatorFactory.hpp"

namespace sgpp {
namespace datadriven {

OperationMultiEvalHPX::OperationMultiEvalHPX(base::Grid& grid,
        base::DataMatrix& dataset,
        sgpp::datadriven::OperationMultipleEvalConfiguration& configuration,
        bool verbose) :
        OperationMultipleEval(grid, dataset), configuration(configuration), dim(
                grid.getDimension()), verbose(verbose), duration(-1.0) {
    // create the kernel specific data structures for the current grid
    this->prepare();
}

OperationMultiEvalHPX::~OperationMultiEvalHPX() {
}

void OperationMultiEvalHPX::mult(sgpp::base::DataVector& alpha,
        sgpp::base::DataVector& result) {
    this->myTimer.start();
    std::cout << "in dummy" << std::endl;

    // create appropriate node level multi eval implementation
    std::unique_ptr<sgpp::base::OperationMultipleEval> nodeMultiEval;
    if (configuration.getType() == OperationMultipleEvalType::STREAMING
            && configuration.getSubType()
                    == OperationMultipleEvalSubType::DEFAULT) {
        nodeMultiEval =
                std::make_unique<datadriven::OperationMultiEvalStreaming>(grid,
                        dataset);
    } else if (configuration.getType() == OperationMultipleEvalType::STREAMING
            && configuration.getSubType()
                    == OperationMultipleEvalSubType::OCLMP) {
        nodeMultiEval = std::unique_ptr<OperationMultipleEval>(
                createStreamingOCLMultiPlatformConfigured(grid, dataset,
                        configuration));
    } else {
        throw base::not_implemented_exception();
    }
    nodeMultiEval->mult(alpha, result);
    this->duration = this->myTimer.stop();
}

void OperationMultiEvalHPX::multTranspose(sgpp::base::DataVector& source,
        sgpp::base::DataVector& result) {
    this->myTimer.start();
    this->duration = this->myTimer.stop();
    throw base::not_implemented_exception();
}

double OperationMultiEvalHPX::getDuration() {
    return this->duration;
}

void OperationMultiEvalHPX::prepare() {
}
}  // namespace datadriven
}  // namespace sgpp
