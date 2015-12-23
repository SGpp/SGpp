// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include "Node.hpp"

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

class OperationOcttreeHistogramRegression {
    base::DataMatrix &dataset;
    base::DataVector &values;
    size_t dims;
public:
    OperationOcttreeHistogramRegression(base::DataMatrix& dataset, base::DataVector &values) :
            dataset(dataset), values(values), dims(dataset.getNcols()) {

    }

    std::unique_ptr<HistogramTree::Node> hierarchize(float_t targetMSE, size_t targetMaxLevel) {

        std::vector<float_t> xRoot(dims);
        for (size_t d = 0; d < dims; d++) {
            xRoot[d] = 0.5;
        }

        std::vector<float_t> hRoot(dims);
        for (size_t d = 0; d < dims; d++) {
            hRoot[d] = 0.5;
        }

        std::vector<size_t> rootSupport(dataset.getNrows());
        for (size_t i = 0; i < dataset.getNrows(); i++) {
            rootSupport[i] = i;
        }

        std::unique_ptr<HistogramTree::Node> root = std::make_unique<HistogramTree::Node>(xRoot, hRoot, rootSupport, dataset, values);

        root->hierarchize(targetMSE, targetMaxLevel);

        return std::move(root);
    }

};

}
}
