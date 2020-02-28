// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <string>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

namespace sgpp {
namespace datadriven {
namespace StreamingModOCLFastMultiPlatform {

class Configuration {
 private:
  Configuration() = default;

 public:
  static const std::string &getKernelName() {
    static std::string kernelName = "StreamingModOCLFastMultiPlatform";
    return kernelName;
  }

  static void augmentDefaultParameters(sgpp::base::OCLOperationConfiguration &parameters) {
    // setup verbose variable for the operation
    if (parameters.contains("VERBOSE") == false) {
      parameters.addIDAttr("VERBOSE", false);
    }

    for (std::string &platformName : parameters["PLATFORMS"].keys()) {
      json::Node &platformNode = parameters["PLATFORMS"][platformName];
      for (std::string &deviceName : platformNode["DEVICES"].keys()) {
        json::Node &deviceNode = platformNode["DEVICES"][deviceName];

        const std::string &kernelName =
            StreamingModOCLFastMultiPlatform::Configuration::getKernelName();

        json::Node &kernelNode = deviceNode["KERNELS"].contains(kernelName)
                                     ? deviceNode["KERNELS"][kernelName]
                                     : deviceNode["KERNELS"].addDictAttr(kernelName);

        // set verbosity for the individual kernels
        if (kernelNode.contains("VERBOSE") == false) {
          kernelNode.addIDAttr("VERBOSE", false);
        }

        if (kernelNode.contains("REUSE_SOURCE") == false) {
          kernelNode.addIDAttr("REUSE_SOURCE", false);
        }

        if (kernelNode.contains("WRITE_SOURCE") == false) {
          kernelNode.addIDAttr("WRITE_SOURCE", false);
        }

        if (kernelNode.contains("LOCAL_SIZE") == false) {
          kernelNode.addIDAttr("LOCAL_SIZE", UINT64_C(128));
        }

        if (kernelNode.contains("KERNEL_SCHEDULE_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_SCHEDULE_SIZE", UINT64_C(102400));
        }

        if (kernelNode.contains("INTERNAL_PRECISION") == false) {
          kernelNode.addTextAttr("INTERNAL_PRECISION", "double");
        }

        if (kernelNode.contains("KERNEL_USE_LOCAL_MEMORY") == false) {
          kernelNode.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
        }

        if (kernelNode.contains("KERNEL_STORE_DATA") == false) {
          kernelNode.addTextAttr("KERNEL_STORE_DATA", "array");
        }

        if (kernelNode.contains("KERNEL_MAX_DIM_UNROLL") == false) {
          kernelNode.addIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(1));
        }

        if (kernelNode.contains("KERNEL_DATA_BLOCK_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_DATA_BLOCK_SIZE", UINT64_C(1));
        }

        if (kernelNode.contains("KERNEL_TRANS_DATA_BLOCK_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_TRANS_DATA_BLOCK_SIZE", UINT64_C(1));
        }

        if (kernelNode.contains("KERNEL_TRANS_GRID_BLOCK_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(1));
        }

        if (kernelNode.contains("KERNEL_PREFETCH_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_PREFETCH_SIZE", UINT64_C(64));
        }
      }
    }
  }
};
}  // namespace StreamingModOCLFastMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
