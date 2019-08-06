// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/datadriven/operation/hash/OperationPiecewiseConstantRegression/Node.hpp>
#include <sgpp/globaldef.hpp>

#include <omp.h>

#include <vector>

namespace sgpp {
namespace datadriven {

class OperationPiecewiseConstantRegression {
  base::DataMatrix& dataset;
  base::DataVector& values;
  size_t dims;

 public:
  OperationPiecewiseConstantRegression(base::DataMatrix& dataset, base::DataVector& values)
      : dataset(dataset), values(values), dims(dataset.getNcols()) {}

  std::unique_ptr<PiecewiseConstantRegression::Node> hierarchize(double targetMSE,
                                                                 size_t targetMaxLevel) {
    std::vector<double> xRoot(dims);

    for (size_t d = 0; d < dims; d++) {
      xRoot[d] = 0.5;
    }

    std::vector<double> hRoot(dims);

    for (size_t d = 0; d < dims; d++) {
      hRoot[d] = 0.5;
    }

    std::vector<size_t> rootSupport(dataset.getNrows());

    for (size_t i = 0; i < dataset.getNrows(); i++) {
      rootSupport[i] = i;
    }

    std::unique_ptr<PiecewiseConstantRegression::Node> root =
        std::make_unique<PiecewiseConstantRegression::Node>(xRoot, hRoot, rootSupport, dataset,
                                                            values);

    root->hierarchize(targetMSE, targetMaxLevel);

    std::cout << "total node count: " << (root->getChildCount() + 1) << std::endl;
    std::cout << "hierarchization max level: " << root->getHierarchizationMaxLevel() << std::endl;

    return root;
  }
};
}  // namespace datadriven
}  // namespace sgpp
