/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include <iostream>

#include <sgpp/base/exception/operation_exception.hpp>

#include "StreamingModOCLFastMultiPlatformKernelSourceBuilder.hpp"

namespace SGPP {
namespace datadriven {

StreamingModOCLFastMultiPlatformKernelSourceBuilder::StreamingModOCLFastMultiPlatformKernelSourceBuilder(
        std::shared_ptr<base::OCLOperationConfiguration> parameters, size_t dims, json::Node &firstDeviceConfig) :
        parameters(parameters), dims(dims), indent("    "), indent2("        "), indent3("            "), indent4(
                "                "), firstDeviceConfig(firstDeviceConfig) {
    localWorkgroupSize = firstDeviceConfig["LOCAL_SIZE"].getUInt();
    useLocalMemory = firstDeviceConfig["KERNEL_USE_LOCAL_MEMORY"].getBool();
    dataBlockSize = firstDeviceConfig["KERNEL_DATA_BLOCKING_SIZE"].getUInt();
    maxDimUnroll = firstDeviceConfig["KERNEL_MAX_DIM_UNROLL"].getUInt();
}

std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::getData(std::string dim, size_t dataBlockingIndex) {
    std::stringstream output;
    if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("array") == 0) {
        output << "data_" << dataBlockingIndex << "[" << dim << "]";
    } else if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("register") == 0) {
        output << "data_" << dataBlockingIndex << "_" << dim;
    } else if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
        output << "ptrData[(" << dataBlockSize << " * globalIdx) + (resultSize * " << dim << ") + " << dataBlockingIndex
                << "]";
    } else {
        throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
}

std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::getData(size_t dim, size_t dataBlockingIndex) {
    std::stringstream dimStringStream;
    dimStringStream << dim;
    std::string dimString = dimStringStream.str();
    return this->getData(dimString, dataBlockingIndex);
}

std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::asString() {
    if (firstDeviceConfig["INTERNAL_PRECISION"].get() == "float") {
        return "float";
    } else {
        return "double";
    }
}
std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::constSuffix() {
    if (firstDeviceConfig["INTERNAL_PRECISION"].get() == "float") {
        return "f";
    } else {
        return "";
    }
}
std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::intAsString() {
    if (firstDeviceConfig["INTERNAL_PRECISION"].get() == "float") {
        return "uint";
    } else {
        return "ulong";
    }
}

std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::unrolledBasisFunctionEvalulation(size_t dims,
        size_t startDim, size_t endDim, std::string unrollVariable) {
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
        if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("register") == 0) {
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

        output << indent3 << "dimLevelIndex = " << "(m * " << dims << ") + " << pointerAccess << ";" << std::endl;

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

std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::getDataTrans(std::string dim,
        size_t dataBlockingIndex) {
    std::stringstream output;
    output << "ptrData[dimDataIndex + " << (localWorkgroupSize * dataBlockingIndex) << "]";
    return output.str();
}

std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::getLevelTrans(std::string dim,
        size_t gridBlockingIndex) {
    std::stringstream output;
    if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("array") == 0) {
        output << "level_" << gridBlockingIndex << "[" << dim << "]";
    } else if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("register") == 0) {
        output << "level_" << gridBlockingIndex << "_" << dim;
    } else if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
        output << "ptrLevel[dimLevelIndex]";
    } else {
        throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
}

std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::getIndexTrans(std::string dim,
        size_t gridBlockingIndex) {
    std::stringstream output;
    if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("array") == 0) {
        output << "index_" << gridBlockingIndex << "[" << dim << "]";
    } else if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("register") == 0) {
        output << "index_" << gridBlockingIndex << "_" << dim;
    } else if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
        output << "ptrIndex[dimLevelIndex]";
    } else {
        throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
}

std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::unrolledBasisFunctionEvalulationTrans(size_t dims,
        size_t startDim, size_t endDim, std::string unrollVariable, size_t gridBlockIndex) {

    size_t transDataBlockSize = firstDeviceConfig["KERNEL_TRANS_DATA_BLOCK_SIZE"].getUInt();
    size_t transGridBlockSize = firstDeviceConfig["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt();

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
        if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("register") == 0) {
            std::stringstream stream;
            stream << (d);
            dString = stream.str();
        } else {
            dString = pointerAccess;
        }

        output << indent3 << "dimDataIndex = " << "(" << pointerAccess << " * sourceSize) + k;" << std::endl;

        if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
            output << indent3 << "dimLevelIndex = " << "((" << transGridBlockSize << " * groupIdx + " << gridBlockIndex
                    << ") * " << dims << ") +" << pointerAccess << ";" << std::endl;
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

std::string StreamingModOCLFastMultiPlatformKernelSourceBuilder::reuseSource(std::string fileName) {
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

void StreamingModOCLFastMultiPlatformKernelSourceBuilder::writeSource(std::string fileName, std::string source) {
    //update file with kernel (for debugging)
    std::ofstream sourceFile;
    sourceFile.open(fileName);
    sourceFile << source;
    sourceFile.close();
}

}
}
