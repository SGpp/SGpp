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

StreamingOCLKernelSourceBuilder::StreamingOCLKernelSourceBuilder(base::OCLConfigurationParameters parameters,
        size_t dims) :
        parameters(parameters), dims(dims), indent("    "), indent2("        "), indent3("            "), indent4(
                "                ") {
    localWorkgroupSize = parameters.getAsUnsigned("LOCAL_SIZE");
    useLocalMemory = this->parameters.getAsBoolean("KERNEL_USE_LOCAL_MEMORY");
    dataBlockSize = parameters.getAsUnsigned("KERNEL_DATA_BLOCKING_SIZE");
    maxDimUnroll = this->parameters.getAsUnsigned("KERNEL_MAX_DIM_UNROLL");
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

std::string StreamingOCLKernelSourceBuilder::getData(std::string dim, size_t dataBlockingIndex) {
    std::stringstream output;
    if (parameters["KERNEL_STORE_DATA"].compare("array") == 0) {
        output << "data_" << dataBlockingIndex << "[" << dim << "]";
    } else if (parameters["KERNEL_STORE_DATA"].compare("register") == 0) {
        output << "data_" << dataBlockingIndex << "_" << dim;
    } else if (parameters["KERNEL_STORE_DATA"].compare("pointer") == 0) {
        output << "ptrData[(" << dataBlockSize << " * globalIdx) + (resultSize * " << dim << ") + " << dataBlockingIndex
                << "]";
    } else {
        throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
}

std::string StreamingOCLKernelSourceBuilder::getData(size_t dim, size_t dataBlockingIndex) {
    std::stringstream dimStringStream;
    dimStringStream << dim;
    std::string dimString = dimStringStream.str();
    return this->getData(dimString, dataBlockingIndex);
}

std::string StreamingOCLKernelSourceBuilder::unrolledBasisFunctionEvalulation(size_t dims, size_t startDim,
        size_t endDim, std::string unrollVariable) {
    std::stringstream output;

    for (size_t d = startDim; d < endDim; d++) {

        std::stringstream dimElement;
        dimElement << "(";
        if (!unrollVariable.compare("") == 0) {
            dimElement << unrollVariable << " + ";
        }
        dimElement << d;
        dimElement << ")";
        std::string pointerAccess = dimElement.str();

        std::string dString;
        if (parameters["KERNEL_STORE_DATA"].compare("register") == 0) {
            std::stringstream stream;
            stream << (d);
            dString = stream.str();
        } else {
            dString = pointerAccess;
        }

        std::stringstream levelAccessStream;
        std::stringstream indexAccessStream;
        if (useLocalMemory) {
            levelAccessStream << "locLevel[dimLevelIndex]";
            indexAccessStream << "locIndex[dimLevelIndex]";
        } else {
            levelAccessStream << "ptrLevel[dimLevelIndex]";
            indexAccessStream << "ptrIndex[dimLevelIndex]";
        }
        std::string levelAccess = levelAccessStream.str();
        std::string indexAccess = indexAccessStream.str();

        output << indent3 << "dimLevelIndex = " << "(k * " << dims << ") + " << pointerAccess << ";" << std::endl;

        for (size_t i = 0; i < dataBlockSize; i++) {
            output << indent3 << "curSupport_" << i << " *= fmax(1.0" << this->constSuffix() << " - fabs((";
            output << levelAccess << " * " << getData(dString, i) << ") - " << indexAccess << "), 0.0"
                    << this->constSuffix() << ");" << std::endl;
        }

//        output << "     eval = (" << levelAccess << " * " << getData(dString, 0) << ");" << std::endl;
//        output << "     index_calc = eval - " << indexAccess << ";" << std::endl;
//        output << "     abs = fabs(index_calc);" << std::endl;
//        output << "     last = 1.0" << this->constSuffix() << " - abs;" << std::endl;
//        output << "     localSupport = fmax(last, 0.0" << this->constSuffix() << ");" << std::endl;
//        output << "     curSupport *= localSupport;" << std::endl << std::endl;

    }
    return output.str();
}

}
}

