/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#pragma once

#include <fstream>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>

namespace SGPP {
namespace datadriven {
namespace StreamingOCLMultiPlatform {

template<typename real_type>
class SourceBuilderMultTranspose: public base::KernelSourceBuilderBase<real_type> {
private:

    std::shared_ptr<base::OCLDevice> device;

    json::Node &kernelConfiguration;

    size_t dims;

    size_t localWorkgroupSize;
    bool useLocalMemory;
    size_t dataBlockSize;
    size_t transGridBlockSize;
    uint64_t maxDimUnroll;

    std::string getLevel(std::string dim, size_t gridBlockingIndex) {
        std::stringstream output;
        if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
            output << "level_" << gridBlockingIndex << "[" << dim << "]";
        } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
            output << "level_" << gridBlockingIndex << "_" << dim;
        } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
            output << "ptrLevel[dimLevelIndex]";
        } else {
            throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
        }
        return output.str();
    }

    std::string getIndex(std::string dim, size_t gridBlockingIndex) {
        std::stringstream output;
        if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
            output << "index_" << gridBlockingIndex << "[" << dim << "]";
        } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
            output << "index_" << gridBlockingIndex << "_" << dim;
        } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
            output << "ptrIndex[dimLevelIndex]";
        } else {
            throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
        }
        return output.str();
    }

    std::string getData(std::string dim, size_t dataBlockingIndex) {
        std::stringstream output;
        if (kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool()) {
            // (locData[(d * 128)+k])
            output << "locData[(" << dim << " * " << localWorkgroupSize << ")+k]";
        } else {
            output << "ptrData[(" << dim << "* sourceSize) + k]";
            //        output << "ptrData[dimDataIndex + " << (localWorkgroupSize * dataBlockingIndex) << "]";
        }
        return output.str();
    }

    std::string unrolledBasisFunctionEvalulation(size_t dims, size_t startDim, size_t endDim,
            std::string unrollVariable) {
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
            if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
                std::stringstream stream;
                stream << (d);
                dString = stream.str();
            } else {
                dString = pointerAccess;
            }

            //        output << indent3 << "dimLevelIndex = " << "(k * " << dims << ") + " << pointerAccess << ";" << std::endl;

            for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
                output << this->indent3 << "curSupport_" << gridIndex << " *= fmax(1.0" << this->constSuffix()
                        << " - fabs((";
                output << getLevel(dString, gridIndex) << " * " << getData(dString, 0) << ") - "
                        << getIndex(dString, gridIndex) << "), 0.0" << this->constSuffix() << ");" << std::endl;
            }
        }
        return output.str();
    }

public:

    SourceBuilderMultTranspose(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration, size_t dims) :
            device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
        localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
        useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
        dataBlockSize = kernelConfiguration["KERNEL_DATA_BLOCKING_SIZE"].getUInt();
        transGridBlockSize = kernelConfiguration["KERNEL_TRANS_GRID_BLOCKING_SIZE"].getUInt();
        maxDimUnroll = kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt();
    }

    std::string generateSource() {

        if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
            return this->reuseSource("StreamingOCLMultiPlatform_multTrans.cl");
        }

        size_t localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
        bool useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
        uint64_t maxDimUnroll = kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt();

        std::stringstream sourceStream;

        if (this->floatType().compare("double") == 0) {
            sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
        }

        sourceStream << "__kernel" << std::endl;

        sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))" << std::endl;
        sourceStream << "void multTransOCL(__global const " << this->floatType() << "* ptrLevel," << std::endl;
        sourceStream << "                  __global const " << this->floatType() << "* ptrIndex," << std::endl;
        sourceStream << "                  __global const " << this->floatType() << "* ptrData," << std::endl;
        sourceStream << "                  __global const " << this->floatType() << "* ptrSource," << std::endl;
        sourceStream << "                  __global       " << this->floatType() << "* ptrResult," << std::endl;
        sourceStream << "                  uint sourceSize," << std::endl;
        sourceStream << "                  uint start_data," << std::endl;
        sourceStream << "                  uint end_data) {" << std::endl;
        sourceStream << this->indent << "int globalIdx = get_global_id(0);" << std::endl;
        sourceStream << this->indent << "int localIdx = get_local_id(0);" << std::endl;
        sourceStream << std::endl;

        for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
            sourceStream << this->indent << this->floatType() << " myResult_" << gridIndex << " = 0.0;" << std::endl;
        }
        sourceStream << std::endl;

        if (useLocalMemory) {
            sourceStream << this->indent << "__local " << this->floatType() << " locData[" << dims * localWorkgroupSize
                    << "];" << std::endl;
            sourceStream << this->indent << "__local " << this->floatType() << " locSource[" << localWorkgroupSize
                    << "];" << std::endl << std::endl;
        }

        // create a register storage for the level and index of the grid points of the work item
        if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
            for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
                sourceStream << this->indent << this->floatType() << " level_" << gridIndex << "[" << dims << "];"
                        << std::endl;
                sourceStream << this->indent << this->floatType() << " index_" << gridIndex << "[" << dims << "];"
                        << std::endl;

                for (size_t d = 0; d < dims; d++) {
                    sourceStream << this->indent << "level_" << gridIndex << "[" << d << "] = ptrLevel[(("
                            << transGridBlockSize << " * globalIdx + " << gridIndex << ") * " << dims << ")+" << d
                            << "];" << std::endl;
                    sourceStream << this->indent << "index_" << gridIndex << "[" << d << "] = ptrIndex[(("
                            << transGridBlockSize << " * globalIdx + " << gridIndex << ") * " << dims << ")+" << d
                            << "];" << std::endl;

                }
                sourceStream << std::endl;
            }
            sourceStream << std::endl;
        } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
            for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {

                for (size_t d = 0; d < dims; d++) {
                    sourceStream << this->indent << this->floatType() << " level_" << gridIndex << "_" << d
                            << " = ptrLevel[((" << transGridBlockSize << " * globalIdx + " << gridIndex << ") * "
                            << dims << ")+" << d << "];" << std::endl;
                    sourceStream << this->indent << this->floatType() << " index_" << gridIndex << "_" << d
                            << " = ptrIndex[((" << transGridBlockSize << " * globalIdx + " << gridIndex << ") * "
                            << dims << ")+" << d << "];" << std::endl;
                }
                sourceStream << std::endl;
            }
            sourceStream << std::endl;
        }

        sourceStream << this->indent << "// Iterate over all data points" << std::endl;

        if (useLocalMemory) {
            sourceStream << this->indent << "for(int i = start_data; i < end_data; i+=" << localWorkgroupSize << ") {"
                    << std::endl;

            if (dims > maxDimUnroll) {
                sourceStream << this->indent2 << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
                sourceStream << this->indent3 << "locData[(d * " << localWorkgroupSize
                        << ")+(localIdx)] = ptrData[(d * sourceSize) + (localIdx + i)];" << std::endl;
                sourceStream << this->indent2 << "}" << std::endl;
            } else {
                for (size_t d = 0; d < dims; d++) {
                    sourceStream << this->indent2 << "locData[(" << d << "*" << localWorkgroupSize
                            << ")+(localIdx)] = ptrData[(" << d << "*sourceSize)+(localIdx+i)];" << std::endl;
                }
            }

            sourceStream << this->indent2 << "locSource[localIdx] = ptrSource[i+localIdx];" << std::endl;
            sourceStream << this->indent2 << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl << std::endl;
            sourceStream << this->indent2 << "for(int k = 0; k < " << localWorkgroupSize << "; k++) {" << std::endl;

            for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
                sourceStream << this->indent3 << this->floatType() << " curSupport_" << gridIndex << " = locSource[k];"
                        << std::endl;
            }
        } else {
            sourceStream << this->indent << "for(int k = start_data; k < end_data; k++) {" << std::endl;

            for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
                sourceStream << this->indent2 << this->floatType() << " curSupport_" << gridIndex << " = ptrSource[k];"
                        << std::endl;
            }
        }
        sourceStream << std::endl;

        if (dims > maxDimUnroll) {
            sourceStream << this->indent2 << "for (size_t unrollDim = 0; unrollDim < "
                    << ((dims / maxDimUnroll) * maxDimUnroll) << "; unrollDim += " << maxDimUnroll << ") {"
                    << std::endl;

            sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, std::min(maxDimUnroll, dims), "unrollDim");
            sourceStream << this->indent2 << "}" << std::endl;

            if (dims % maxDimUnroll != 0) {
                sourceStream
                        << this->unrolledBasisFunctionEvalulation(dims, (dims / maxDimUnroll) * maxDimUnroll, dims, "");
            }

        } else {
            sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, dims, "");
        }

        sourceStream << std::endl;

        for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
            sourceStream << this->indent2 << "myResult_" << gridIndex << " += curSupport_" << gridIndex << ";"
                    << std::endl;
        }
        sourceStream << this->indent << "}" << std::endl << std::endl;

        if (useLocalMemory) {
            sourceStream << this->indent2 << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
            sourceStream << this->indent << "}" << std::endl;
        }

        for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {

            sourceStream << this->indent << "ptrResult[(" << transGridBlockSize << " * globalIdx) + " << gridIndex
                    << "] = myResult_" << gridIndex << ";" << std::endl;
        }
        sourceStream << "}" << std::endl;

        if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
            this->writeSource("StreamingOCLMultiPlatform_multTrans.cl", sourceStream.str());
        }

        return sourceStream.str();
    }

}
;

}
}
}
