/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include <iostream>

#include <sgpp/base/exception/operation_exception.hpp>

#include "StreamingModOCLFastKernelSourceBuilder.hpp"

namespace SGPP {
namespace datadriven {

StreamingModOCLFastKernelSourceBuilder::StreamingModOCLFastKernelSourceBuilder(base::ConfigurationParameters parameters,
        size_t dims) :
        parameters(parameters), dims(dims), indent("    "), indent2("        "), indent3("            "), indent4(
                "                ") {
    localWorkgroupSize = parameters.getAsUnsigned("LOCAL_SIZE");
    useLocalMemory = parameters.getAsBoolean("KERNEL_USE_LOCAL_MEMORY");
    dataBlockSize = parameters.getAsUnsigned("KERNEL_DATA_BLOCKING_SIZE");
    maxDimUnroll = parameters.getAsUnsigned("KERNEL_MAX_DIM_UNROLL");
}

std::string StreamingModOCLFastKernelSourceBuilder::getData(std::string dim, size_t dataBlockingIndex) {
    std::stringstream output;
    if (parameters["KERNEL_STORE_DATA"].compare("array") == 0) {
        output << "data_" << dataBlockingIndex << "[" << dim << "]";
    } else if (parameters["KERNEL_STORE_DATA"].compare("pointer") == 0) {
        output << "ptrData[(" << dataBlockSize << " * globalIdx) + (resultSize * " << dim << ") + " << dataBlockingIndex
                << "]";
    } else {
        throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
}

std::string StreamingModOCLFastKernelSourceBuilder::getData(size_t dim, size_t dataBlockingIndex) {
    std::stringstream dimStringStream;
    dimStringStream << dim;
    std::string dimString = dimStringStream.str();
    return this->getData(dimString, dataBlockingIndex);
}

std::string StreamingModOCLFastKernelSourceBuilder::asString() {
    if (parameters["INTERNAL_PRECISION"] == "float") {
        return "float";
    } else {
        return "double";
    }
}
std::string StreamingModOCLFastKernelSourceBuilder::constSuffix() {
    if (parameters["INTERNAL_PRECISION"] == "float") {
        return "f";
    } else {
        return "";
    }
}
std::string StreamingModOCLFastKernelSourceBuilder::intAsString() {
    if (parameters["INTERNAL_PRECISION"] == "float") {
        return "uint";
    } else {
        return "ulong";
    }
}

std::string StreamingModOCLFastKernelSourceBuilder::unrolledBasisFunctionEvalulation(size_t dims, size_t startDim,
        size_t endDim, std::string unrollVariable) {
    std::stringstream output;

    for (size_t d = startDim; d < endDim; d++) {

        std::stringstream dimElement;
        if (!unrollVariable.compare("") == 0) {
            dimElement << unrollVariable << " + ";
        }
        dimElement << d;
        std::string dString = dimElement.str();

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

        output << indent3 << "dimLevelIndex = " << "(m * " << dims << ") + " << dString << ";" << std::endl;

        output << indent3 << "if (" << levelAccess << " == 2.0" << this->constSuffix() << ") {" << std::endl;
        output << indent3 << "} else if (" << indexAccess << " == 1.0" << this->constSuffix() << ") {" << std::endl;

        for (size_t i = 0; i < dataBlockSize; i++) {
            output << indent4 << "curSupport_" << i << " *= max(2.0" << this->constSuffix() << " - (" << levelAccess
                    << " * " << getData(dString, i) << "), 0.0" << this->constSuffix() << ") ;" << std::endl;
        }

        output << indent3 << "} else if (" << indexAccess << " == " << levelAccess << " - 1.0" << this->constSuffix()
                << ")" << std::endl;
        output << indent3 << "{" << std::endl;

        for (size_t i = 0; i < dataBlockSize; i++) {
            output << indent4 << "curSupport_" << i << " *= max((" << levelAccess << " * " << getData(dString, i)
                    << ") - " << indexAccess << " + 1.0" << this->constSuffix() << ", 0.0" << this->constSuffix()
                    << ");" << std::endl;
        }

        output << indent3 << "} else {" << std::endl;

        for (size_t i = 0; i < dataBlockSize; i++) {
            output << indent4 << "curSupport_" << i << " *= max(1.0" << this->constSuffix() << " - fabs(" << levelAccess
                    << " * " << getData(dString, i) << " - " << indexAccess << "), 0.0" << this->constSuffix() << ");"
                    << std::endl;
        }
        output << indent3 << "}" << std::endl;
    }
    return output.str();
}

std::string StreamingModOCLFastKernelSourceBuilder::getDataTrans(std::string dim, size_t dataBlockingIndex) {
    std::stringstream output;
    output << "ptrData[dimDataIndex + " << (localWorkgroupSize * dataBlockingIndex) << "]";
    return output.str();
}

std::string StreamingModOCLFastKernelSourceBuilder::getLevelTrans(std::string dim, size_t gridBlockingIndex) {
    std::stringstream output;
    if (parameters["KERNEL_STORE_DATA"].compare("array") == 0) {
        output << "level_" << gridBlockingIndex << "[" << dim << "]";
    } else if (parameters["KERNEL_STORE_DATA"].compare("pointer") == 0) {
        output << "ptrLevel[dimLevelIndex]";
    } else {
        throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
}

std::string StreamingModOCLFastKernelSourceBuilder::getIndexTrans(std::string dim, size_t gridBlockingIndex) {
    std::stringstream output;
    if (parameters["KERNEL_STORE_DATA"].compare("array") == 0) {
        output << "index_" << gridBlockingIndex << "[" << dim << "]";
    } else if (parameters["KERNEL_STORE_DATA"].compare("pointer") == 0) {
        output << "ptrIndex[dimLevelIndex]";
    } else {
        throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
}

std::string StreamingModOCLFastKernelSourceBuilder::unrolledBasisFunctionEvalulationTrans(size_t dims, size_t startDim,
        size_t endDim, std::string unrollVariable, size_t gridBlockIndex) {

    size_t transDataBlockSize = this->parameters.getAsUnsigned("KERNEL_TRANS_DATA_BLOCK_SIZE");
    size_t transGridBlockSize = this->parameters.getAsUnsigned("KERNEL_TRANS_GRID_BLOCK_SIZE");

    std::stringstream output;

    for (size_t d = startDim; d < endDim; d++) {

        std::stringstream dimElement;
        dimElement << "(";
        if (!unrollVariable.compare("") == 0) {
            dimElement << unrollVariable << " + ";
        }
        dimElement << d;
        dimElement << ")";
        std::string dString = dimElement.str();

        output << indent3 << "dimDataIndex = " << "(" << dString << " * sourceSize) + k;" << std::endl;

        if (parameters["KERNEL_STORE_DATA"].compare("pointer") == 0) {
            output << indent3 << "dimLevelIndex = " << "((" << transGridBlockSize << " * groupIdx + " << gridBlockIndex
                    << ") * " << dims << ") +" << dString << ";" << std::endl;
        }

        output << indent3 << "if (" << getLevelTrans(dString, gridBlockIndex) << " == 2.0" << this->constSuffix()
                << ") {" << std::endl;

        output << indent3 << "} else if (" << getIndexTrans(dString, gridBlockIndex) << " == 1.0" << this->constSuffix()
                << ") {" << std::endl;

        for (size_t dataIndex = 0; dataIndex < transDataBlockSize; dataIndex++) {
            output << indent4 << "curSupport_" << gridBlockIndex << "_" << dataIndex << " *= max(2.0"
                    << this->constSuffix() << " - (" << getLevelTrans(dString, gridBlockIndex) << " * "
                    << getDataTrans(dString, dataIndex) << "), 0.0" << this->constSuffix() << ") ;" << std::endl;
        }

        output << indent3 << "} else if (" << getIndexTrans(dString, gridBlockIndex) << " == ("
                << getLevelTrans(dString, gridBlockIndex) << " - 1.0" << this->constSuffix() << ")) {" << std::endl;

        for (size_t dataIndex = 0; dataIndex < transDataBlockSize; dataIndex++) {
            output << indent4 << "curSupport_" << gridBlockIndex << "_" << dataIndex << " *= max(("
                    << getLevelTrans(dString, gridBlockIndex) << " * " << getDataTrans(dString, dataIndex) << ") - "
                    << getIndexTrans(dString, gridBlockIndex) << " + 1.0, 0.0);" << std::endl;
        }

        output << indent3 << "} else {" << std::endl;

        for (size_t dataIndex = 0; dataIndex < transDataBlockSize; dataIndex++) {
            output << indent4 << "curSupport_" << gridBlockIndex << "_" << dataIndex << " *= max(1.0"
                    << this->constSuffix() << " - fabs(" << getLevelTrans(dString, gridBlockIndex) << " * "
                    << getDataTrans(dString, dataIndex) << " - " << getIndexTrans(dString, gridBlockIndex) << "), 0.0"
                    << this->constSuffix() << ");" << std::endl;
        }

        output << indent3 << "}" << std::endl;
    }
    return output.str();
}

std::string StreamingModOCLFastKernelSourceBuilder::reuseSource(std::string fileName) {
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

void StreamingModOCLFastKernelSourceBuilder::writeSource(std::string fileName, std::string source) {
    //update file with kernel (for debugging)
    std::ofstream sourceFile;
    sourceFile.open(fileName);
    sourceFile << source;
    sourceFile.close();
}

}
}
