/*
 * StreamingOCLMultiPlatformConfiguration.hpp
 *
 *  Created on: Nov 18, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

class StreamingModOCLFastMultiPlatformConfiguration {
private:
    StreamingModOCLFastMultiPlatformConfiguration() = default;
public:
    static const std::string &getKernelName() {
        static std::string kernelName = "StreamingModOCLFastMultiPlatform";
        return kernelName;
    }

    static void augmentDefaultParameters(SGPP::base::OCLOperationConfiguration &parameters) {

        for (std::string &platformName : parameters["PLATFORMS"].keys()) {
            json::Node &platformNode = parameters["PLATFORMS"][platformName];
            for (std::string &deviceName : platformNode["DEVICES"].keys()) {
                json::Node &deviceNode = platformNode["DEVICES"][deviceName];

                const std::string &kernelName = SGPP::datadriven::StreamingModOCLFastMultiPlatformConfiguration::getKernelName();

                json::Node &kernelNode =
                        deviceNode["KERNELS"].contains(kernelName) ?
                                deviceNode["KERNELS"][kernelName] : deviceNode["KERNELS"].addDictAttr(kernelName);
    //            std::cout << "in kernel augment" << std::endl;
    //            std::cout << "-----------------------------------" << std::endl;
    //            for (std::string &key: kernelNode.keys()) {
    //                std::cout << "key: " << key << " value: " << kernelNode[key].get() << std::endl;
    //            }

                if (kernelNode.contains("LOCAL_SIZE") == false) {
                    kernelNode.addIDAttr("LOCAL_SIZE", 128ul);
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
                    kernelNode.addIDAttr("KERNEL_MAX_DIM_UNROLL", 10ul);
                }

                if (kernelNode.contains("KERNEL_DATA_BLOCKING_SIZE") == false) {
                    kernelNode.addIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
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

}
}
