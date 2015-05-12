/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include "StreamingOCLKernelSourceBuilder.hpp"

#include <fstream>
#include <sstream>

#include "StreamingOCLKernelSourceBuilder.hpp"

#include <sgpp/base/exception/operation_exception.hpp>

namespace SGPP {
namespace datadriven {

StreamingOCLKernelSourceBuilder::StreamingOCLKernelSourceBuilder(base::OCLConfigurationParameters parameters) :
        parameters(parameters) {

}

std::string StreamingOCLKernelSourceBuilder::asString() {
    if (parameters["INTERNAL_PRECISION"] == "float") {
        return "float";
    } else {
        return "double";
    }
}
std::string StreamingOCLKernelSourceBuilder::constSuffix() {
    if (parameters["INTERNAL_PRECISION"] == "float") {
        return "f";
    } else {
        return "";
    }
}
std::string StreamingOCLKernelSourceBuilder::intAsString() {
    if (parameters["INTERNAL_PRECISION"] == "float") {
        return "uint";
    } else {
        return "ulong";
    }
}

std::string StreamingOCLKernelSourceBuilder::reuseSource(std::string fileName) {
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

void StreamingOCLKernelSourceBuilder::writeSource(std::string fileName, std::string source) {
    //update file with kernel (for debugging)
    std::ofstream sourceFile;
    sourceFile.open(fileName);
    sourceFile << source;
    sourceFile.close();
}

}
}

