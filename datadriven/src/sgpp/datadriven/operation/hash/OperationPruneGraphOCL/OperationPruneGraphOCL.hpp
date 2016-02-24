// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>

#include <vector>

namespace SGPP {
namespace datadriven {
namespace DensityOCLMultiPlatform {

class OperationPruneGraphOCL {
 public:
    OperationPruneGraphOCL()  {
    }

    virtual void prune_graph(std::vector<int> &graph) = 0;
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace SGPP
