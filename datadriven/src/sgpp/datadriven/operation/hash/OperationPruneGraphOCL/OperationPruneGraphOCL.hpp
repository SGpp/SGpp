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
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/KernelPruneGraph.hpp>

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// Pure virtual base class for the graph pruning operation
class OperationPruneGraphOCL {
 public:
  OperationPruneGraphOCL()  {
  }


  /// Deletes all nodes and edges within areas of low density which are in the given graph chunk
  virtual void prune_graph(std::vector<int> &graph, size_t startid = 0, size_t chunksize = 0) = 0;
  virtual ~OperationPruneGraphOCL(void) {}
  static void load_default_parameters(base::OCLOperationConfiguration *parameters) {
    if (parameters->contains("INTERNAL_PRECISION") == false) {
      std::cout << "Warning! No internal precision setting detected."
                << " Using double precision from now on!" << std::endl;
      parameters->addIDAttr("INTERNAL_PRECISION", "double");
    }
    if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
      DensityOCLMultiPlatform::KernelPruneGraph<float>::augmentDefaultParameters(*parameters);
    } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
      DensityOCLMultiPlatform::KernelPruneGraph<double>::augmentDefaultParameters(*parameters);
    } else {
      std::stringstream errorString;
      errorString << "Error creating operation\"CreateGraphOCL\": "
                  << " invalid value for parameter \"INTERNAL_PRECISION\"";
      throw base::operation_exception(errorString.str().c_str());
    }
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
