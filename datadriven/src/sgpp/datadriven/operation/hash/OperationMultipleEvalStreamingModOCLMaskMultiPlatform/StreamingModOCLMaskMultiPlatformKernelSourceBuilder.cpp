/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <fstream>

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/base/grid/GridStorage.hpp>

#include "StreamingModOCLMaskMultiPlatformKernelSourceBuilder.hpp"

namespace SGPP {
namespace datadriven {

StreamingModOCLMaskMultiPlatformKernelSourceBuilder::StreamingModOCLMaskMultiPlatformKernelSourceBuilder(
        std::shared_ptr<base::OCLOperationConfiguration> parameters, size_t dims) :
        parameters(parameters), dims(dims), indent("    "), indent2("        "), indent3("            "), indent4(
                "                ") {
    localWorkgroupSize = (*parameters)["LOCAL_SIZE"].getUInt();
    useLocalMemory = (*parameters)["KERNEL_USE_LOCAL_MEMORY"].getBool();
//    maxDimUnroll = (*parameters)["KERNEL_MAX_DIM_UNROLL"].getUInt();
}

std::string StreamingModOCLMaskMultiPlatformKernelSourceBuilder::asString() {
    if ((*parameters)["INTERNAL_PRECISION"].get() == "float") {
        return "float";
    } else {
        return "double";
    }
}

std::string StreamingModOCLMaskMultiPlatformKernelSourceBuilder::constSuffix() {
    if ((*parameters)["INTERNAL_PRECISION"].get() == "float") {
        return "f";
    } else {
        return "";
    }
}

std::string StreamingModOCLMaskMultiPlatformKernelSourceBuilder::intAsString() {
    if ((*parameters)["INTERNAL_PRECISION"].get() == "float") {
        return "uint";
    } else {
        return "ulong";
    }
}

std::string StreamingModOCLMaskMultiPlatformKernelSourceBuilder::reuseSource(std::string fileName) {
    std::stringstream sourceStream;
    std::ifstream file;
    file.open(fileName);

    if (file.is_open()) {
        std::string line;

        while (getline(file, line)) {
            sourceStream << line << std::endl;
        }

        file.close();
    } else {
        throw new base::operation_exception("OCL error: file to reuse not found\n");
    }

    return sourceStream.str();
}

void StreamingModOCLMaskMultiPlatformKernelSourceBuilder::writeSource(std::string fileName, std::string source) {
    //update file with kernel (for debugging)
    std::ofstream sourceFile;
    sourceFile.open(fileName);
    sourceFile << source;
    sourceFile.close();
}

}
}

