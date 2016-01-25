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
class SourceBuilderMult: public base::KernelSourceBuilderBase<real_type> {
private:

    std::shared_ptr<base::OCLDevice> device;

    json::Node &kernelConfiguration;

    size_t dims;

    size_t localWorkgroupSize;
    bool useLocalMemory;
    size_t dataBlockSize;
    size_t transGridBlockSize;
    uint64_t maxDimUnroll;

    std::string getData(std::string dim, size_t dataBlockingIndex) {
        std::stringstream output;
        if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
            output << "data_" << dataBlockingIndex << "[" << dim << "]";
        } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
            output << "data_" << dataBlockingIndex << "_" << dim;
        } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
            output << "ptrData[(" << dataBlockSize << " * globalIdx) + (resultSize * " << dim << ") + "
                    << dataBlockingIndex << "]";
        } else {
            throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
        }
        return output.str();
    }

    std::string getData(size_t dim, size_t dataBlockingIndex) {
        std::stringstream dimStringStream;
        dimStringStream << dim;
        std::string dimString = dimStringStream.str();
        return this->getData(dimString, dataBlockingIndex);
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

            output << this->indent3 << "dimLevelIndex = " << "(k * " << dims << ") + " << pointerAccess << ";"
                    << std::endl;

            for (size_t i = 0; i < dataBlockSize; i++) {
                output << this->indent3 << "curSupport_" << i << " *= fmax(1.0" << this->constSuffix() << " - fabs((";
                output << levelAccess << " * " << getData(dString, i) << ") - " << indexAccess << "), 0.0"
                        << this->constSuffix() << ");" << std::endl;
            }
        }
        return output.str();
    }

public:

    SourceBuilderMult(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration, size_t dims) :
            device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
        localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
        useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
        dataBlockSize = kernelConfiguration["KERNEL_DATA_BLOCKING_SIZE"].getUInt();
        transGridBlockSize = kernelConfiguration["KERNEL_TRANS_GRID_BLOCKING_SIZE"].getUInt();
        maxDimUnroll = kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt();
    }

    std::string generateSource() {

        if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
            return this->reuseSource("StreamingOCLMultiPlatform_mult.cl");
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
        sourceStream << "void multOCL(__global const " << this->floatType() << "* ptrLevel," << std::endl;
        sourceStream << "             __global const " << this->floatType() << "* ptrIndex," << std::endl;
        sourceStream << "             __global const " << this->floatType() << "* ptrData," << std::endl;
        sourceStream << "             __global const " << this->floatType() << "* ptrAlpha," << std::endl;
        sourceStream << "             __global       " << this->floatType() << "* ptrResult," << std::endl;
        sourceStream << "             uint resultSize," << std::endl;
        sourceStream << "             uint start_grid," << std::endl;
        sourceStream << "             uint end_grid) {" << std::endl;
        sourceStream << this->indent << "int globalIdx = get_global_id(0);" << std::endl;
        sourceStream << this->indent << "int localIdx = get_local_id(0);" << std::endl;

        sourceStream << std::endl;

        if (useLocalMemory) {
            sourceStream << this->indent << "__local " << this->floatType() << " locLevel[" << dims * localWorkgroupSize
                    << "];" << std::endl;
            sourceStream << this->indent << "__local " << this->floatType() << " locIndex[" << dims * localWorkgroupSize
                    << "];" << std::endl;
            sourceStream << this->indent << "__local " << this->floatType() << " locAlpha[" << localWorkgroupSize
                    << "];" << std::endl;
            sourceStream << std::endl;
        }

        //    sourceStream << " " << this->floatType() << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl
        //            << std::endl;

        for (size_t i = 0; i < dataBlockSize; i++) {
            sourceStream << this->indent << this->floatType() << " myResult_" << i << " = 0.0;" << std::endl;
        }
        sourceStream << std::endl;

        //caching data in register array, this also requires loading the data into the registers (in contrast using pointers to data directly)
        if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
            for (size_t i = 0; i < dataBlockSize; i++) {
                sourceStream << this->indent << this->floatType() << " data_" << i << "[" << dims << "];" << std::endl;
            }
            sourceStream << std::endl;
            for (size_t i = 0; i < dataBlockSize; i++) {
                for (size_t d = 0; d < dims; d++) {
                    sourceStream << this->indent << getData(d, i) << " = ptrData[" << i << " + (" << dataBlockSize
                            << " * globalIdx) + (resultSize * " << d << ")];" << std::endl;
                }
                sourceStream << std::endl;
            }
        } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
            for (size_t i = 0; i < dataBlockSize; i++) {
                for (size_t d = 0; d < dims; d++) {
                    sourceStream << this->indent << this->floatType() << " " << getData(d, i) << " = ptrData[" << i
                            << " + (" << dataBlockSize << " * globalIdx) + (resultSize * " << d << ")];" << std::endl;
                }
                sourceStream << std::endl;
            }
        }

        sourceStream << this->indent << "size_t dimLevelIndex;" << std::endl;

        sourceStream << std::endl;

        if (useLocalMemory) {
            sourceStream << " // Iterate over all grid points (fast ones, with cache)" << std::endl;
            sourceStream << this->indent << "uint chunkSizeGrid = end_grid - start_grid;" << std::endl;
            sourceStream << this->indent << "uint fastChunkSizeGrid = (chunkSizeGrid / " << localWorkgroupSize << ") * "
                    << localWorkgroupSize << ";" << std::endl;
            sourceStream << this->indent << "for(int j = start_grid; j < start_grid + fastChunkSizeGrid; j+="
                    << localWorkgroupSize << ") {" << std::endl;

            if (dims > maxDimUnroll) {
                sourceStream << this->indent2 << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
                sourceStream << this->indent3 << "locLevel[(localIdx*" << dims << ")+ d] = ptrLevel[((j+localIdx)*"
                        << dims << ")+ d];" << std::endl;
                sourceStream << this->indent3 << "locIndex[(localIdx*" << dims << ")+ d] = ptrIndex[((j+localIdx)*"
                        << dims << ")+ d];" << std::endl;
                sourceStream << this->indent2 << "}" << std::endl;
            } else {
                for (size_t d = 0; d < dims; d++) {
                    sourceStream << this->indent2 << "locLevel[(localIdx*" << dims << ")+" << d
                            << "] = ptrLevel[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
                    sourceStream << this->indent2 << "locIndex[(localIdx*" << dims << ")+" << d
                            << "] = ptrIndex[((j+localIdx)*" << dims << ")+" << d << "];" << std::endl;
                }
            }

            sourceStream << this->indent2 << "locAlpha[localIdx] = ptrAlpha[j+localIdx];" << std::endl;
            sourceStream << this->indent2 << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
            sourceStream << std::endl;
            sourceStream << this->indent2 << "for(int k = 0; k < " << localWorkgroupSize << "; k++) {" << std::endl;
            for (size_t i = 0; i < dataBlockSize; i++) {
                sourceStream << this->indent3 << this->floatType() << " curSupport_" << i << " = locAlpha[k];"
                        << std::endl;
            }
            sourceStream << std::endl;
        } else {
            sourceStream << this->indent << "for(int k = start_grid; k < end_grid; k++) {" << std::endl;
            for (size_t i = 0; i < dataBlockSize; i++) {
                sourceStream << this->indent2 << this->floatType() << " curSupport_" << i << " = ptrAlpha[k];"
                        << std::endl;
            }
            sourceStream << std::endl;
        }

        //    size_t evenDimensions = (dims / dataBlockSize) * dataBlockSize;
        //    size_t remainingDimensions = evenDimensions - remainingDimensions;
        //    sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, dims, "");

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

        for (size_t i = 0; i < dataBlockSize; i++) {
            sourceStream << this->indent2 << "myResult_" << i << " += curSupport_" << i << ";" << std::endl;
        }
        sourceStream << this->indent << "}" << std::endl;
        sourceStream << std::endl;

        if (useLocalMemory) {
            sourceStream << "  barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
            sourceStream << " }" << std::endl;
        }

        sourceStream << std::endl;
        for (size_t i = 0; i < dataBlockSize; i++) {
            sourceStream << this->indent << "ptrResult[(" << dataBlockSize << " * globalIdx) + " << i << "] = myResult_"
                    << i << ";" << std::endl;
        }
        sourceStream << "}" << std::endl;

        if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
            this->writeSource("StreamingOCLMultiPlatform_mult.cl", sourceStream.str());
        }

        return sourceStream.str();
    }

}
;

}
}
}
