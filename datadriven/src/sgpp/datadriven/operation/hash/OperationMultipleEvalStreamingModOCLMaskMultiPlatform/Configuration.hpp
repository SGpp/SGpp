// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
//
//  Created on: Nov 18, 2015
//      Author: pfandedd
//

#pragma once

#include <string>

#include "sgpp/globaldef.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"

namespace SGPP {
namespace datadriven {
namespace StreamingModOCLMaskMultiPlatform {

class Configuration {
 private:
  Configuration() = delete;

 public:
  static const std::string &getKernelName() {
    static std::string kernelName = "StreamingModOCLMaskMultiPlatform";
    return kernelName;
  }

  static void augmentDefaultParameters(SGPP::base::OCLOperationConfiguration &parameters) {
    for (std::string &platformName : parameters["PLATFORMS"].keys()) {
      json::Node &platformNode = parameters["PLATFORMS"][platformName];
      for (std::string &deviceName : platformNode["DEVICES"].keys()) {
        json::Node &deviceNode = platformNode["DEVICES"][deviceName];

        const std::string &kernelName =
            SGPP::datadriven::StreamingModOCLMaskMultiPlatform::Configuration::getKernelName();

        json::Node &kernelNode = deviceNode["KERNELS"].contains(kernelName)
                                     ? deviceNode["KERNELS"][kernelName]
                                     : deviceNode["KERNELS"].addDictAttr(kernelName);
        //            std::cout << "in kernel augment" << std::endl;
        //            std::cout << "-----------------------------------" <<
        //            std::endl;
        //            for (std::string &key: kernelNode.keys()) {
        //                std::cout << "key: " << key << " value: " <<
        //                kernelNode[key].get() << std::endl;
        //            }

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
          kernelNode.addIDAttr("LOCAL_SIZE", 128ul);
        }

        // TODO(pfandedd): if schedule size is not a multiple of LOCAL_SIZE, the kernel crashes with
        // a cryptic error messag, improve
        if (kernelNode.contains("KERNEL_SCHEDULE_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_SCHEDULE_SIZE", 12800ul);
        }

        // TODO(pfandedd): assumed on top level, make consistent
        //                if (kernelNode.contains("INTERNAL_PRECISION") ==
        //                false) {
        //                    kernelNode.addTextAttr("INTERNAL_PRECISION",
        //                    "double");
        //                }

        if (kernelNode.contains("KERNEL_USE_LOCAL_MEMORY") == false) {
          kernelNode.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
        }

        // TODO(pfandedd): implement these parameters in the kernel
        if (kernelNode.contains("KERNEL_STORE_DATA") == false) {
          kernelNode.addTextAttr("KERNEL_STORE_DATA", "array");
        }

        if (kernelNode.contains("KERNEL_MAX_DIM_UNROLL") == false) {
          kernelNode.addIDAttr("KERNEL_MAX_DIM_UNROLL", 10ul);
        }

        if (kernelNode.contains("KERNEL_DATA_BLOCK_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_DATA_BLOCK_SIZE", 1ul);
        }

        if (kernelNode.contains("KERNEL_TRANS_DATA_BLOCK_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_TRANS_DATA_BLOCK_SIZE", 1ul);
        }

        if (kernelNode.contains("KERNEL_TRANS_GRID_BLOCK_SIZE") == false) {
          kernelNode.addIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", 1ul);
        }
      }
    }
  }
};
}  // namespace StreamingModOCLMaskMultiPlatform
}  // namespace datadriven
}  // namespace SGPP
