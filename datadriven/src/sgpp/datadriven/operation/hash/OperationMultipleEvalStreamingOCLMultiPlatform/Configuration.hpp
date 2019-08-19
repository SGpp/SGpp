// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <string>

#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {
namespace StreamingOCLMultiPlatform {

/**
 * Configuration class for the operation OperationMultiEvalStreamingOCLMultiPlatform.
 * This class provides a name and a default configuration for the operation.
 *
 * @see OperationMultiEvalStreamingOCLMultiPlatform
 */
class Configuration {
 private:
  Configuration() = default;

 public:
  /**
   * @return The name of the kernel
   */
  static const std::string &getKernelName() {
    static std::string kernelName = "StreamingOCLMultiPlatform";
    return kernelName;
  }

  /**
   * Takes a configuration and checks whether parameters are missing.
   * Missing parameters are added and initialized with their default values.
   * This function only configures devices that are specified in the configuration.
   *
   * @param parameters Set of parameters to augment with default parameters
   */
  static void augmentDefaultParameters(sgpp::base::OCLOperationConfiguration &parameters) {
    for (std::string &platformName : parameters["PLATFORMS"].keys()) {
      json::Node &platformNode = parameters["PLATFORMS"][platformName];
      for (std::string &deviceName : platformNode["DEVICES"].keys()) {
        json::Node &deviceNode = platformNode["DEVICES"][deviceName];

        const std::string &kernelName = Configuration::getKernelName();

        json::Node &kernelNode = deviceNode["KERNELS"].contains(kernelName)
                                     ? deviceNode["KERNELS"][kernelName]
                                     : deviceNode["KERNELS"].addDictAttr(kernelName);

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

        if (kernelNode.contains("KERNEL_USE_LOCAL_MEMORY") == false) {
          kernelNode.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
        }

        if (kernelNode.contains("KERNEL_STORE_DATA") == false) {
          kernelNode.addTextAttr("KERNEL_STORE_DATA", "array");
        }

        if (kernelNode.contains("KERNEL_MAX_DIM_UNROLL") == false) {
          kernelNode.addIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
        }

        if (kernelNode.contains("KERNEL_DATA_BLOCK_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_DATA_BLOCK_SIZE", UINT64_C(1));
        }

        if (kernelNode.contains("KERNEL_TRANS_GRID_BLOCK_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(1));
        }

        if (kernelNode.contains("KERNEL_SCHEDULE_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_SCHEDULE_SIZE", UINT64_C(102400));
        }

        if (kernelNode.contains("KERNEL_PREFETCH_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_PREFETCH_SIZE", UINT64_C(64));
        }

        if (kernelNode.contains("KERNEL_TRANS_PREFETCH_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_TRANS_PREFETCH_SIZE", UINT64_C(64));
        }
      }
    }
  }
};
}  // namespace StreamingOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
