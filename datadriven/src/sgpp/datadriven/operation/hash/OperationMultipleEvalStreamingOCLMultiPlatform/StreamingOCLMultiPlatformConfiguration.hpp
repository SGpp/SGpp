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

class StreamingOCLMultiPlatformConfiguration {
private:
    StreamingOCLMultiPlatformConfiguration() = default;
public:
    static const std::string &getKernelName() {
        static std::string kernelName = "StreamingOCLMultiPlatform";
        return kernelName;
    }

};

}
}
