// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/exception/operation_exception.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp"

namespace sgpp {
namespace datadriven {

/**
 * This class is a HPX wrapper for other MultiEval-operations that uses a very simple master-slave
 * distributed processing.
 */
class OperationMultiEvalHPX: public sgpp::base::OperationMultipleEval {
protected:
    OperationMultipleEvalType nodeImplType;
    OperationMultipleEvalSubType nodeImplSubType;

    size_t dim;

    bool verbose;

    double duration;

public:
    OperationMultiEvalHPX(base::Grid& grid, base::DataMatrix& dataset,
            OperationMultipleEvalType type, OperationMultipleEvalSubType,
            bool verbose = false);

    ~OperationMultiEvalHPX();

    void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

    void multTranspose(sgpp::base::DataVector& source,
            sgpp::base::DataVector& result);

    void prepare();

    double getDuration();
};

}  // namespace datadriven
}  // namespace sgpp
