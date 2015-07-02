/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include <sstream>

#include <sgpp/base/exception/operation_exception.hpp>
#include "StreamingOCLMPKernelSourceBuilder.hpp"

namespace SGPP {
namespace datadriven {

std::string StreamingOCLMPKernelSourceBuilder::generateSourceMult() {

    if (parameters.getAsBoolean("REUSE_SOURCE")) {
        return this->reuseSource("StreamingOCLMP_mult.cl");
    }

    size_t localWorkgroupSize = parameters.getAsUnsigned("LOCAL_SIZE");
    bool useLocalMemory = this->parameters.getAsBoolean("KERNEL_USE_LOCAL_MEMORY");
    uint64_t maxDimUnroll = this->parameters.getAsUnsigned("KERNEL_MAX_DIM_UNROLL");

    std::stringstream sourceStream;

    if (this->asString() == "double") {
        sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
    }

    sourceStream << "__kernel" << std::endl;
    sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))" << std::endl;
    sourceStream << "void multOCL(__global const " << this->asString() << "* ptrLevel," << std::endl;
    sourceStream << "           __global const " << this->asString() << "* ptrIndex," << std::endl;
    sourceStream << "           __global const " << this->asString() << "* ptrData," << std::endl;
    sourceStream << "           __global const " << this->asString() << "* ptrAlpha," << std::endl;
    sourceStream << "           __global       " << this->asString() << "* ptrResult," << std::endl;
    sourceStream << "           uint resultSize," << std::endl;
    sourceStream << "           uint start_grid," << std::endl;
    sourceStream << "           uint end_grid) " << std::endl;
    sourceStream << "{" << std::endl;
    sourceStream << " int globalIdx = get_global_id(0);" << std::endl;
    sourceStream << " int localIdx = get_local_id(0);" << std::endl;

    sourceStream << std::endl;

    if (useLocalMemory) {
        sourceStream << " __local " << this->asString() << " locLevel[" << dims * localWorkgroupSize << "];"
                << std::endl;
        sourceStream << " __local " << this->asString() << " locIndex[" << dims * localWorkgroupSize << "];"
                << std::endl;
        sourceStream << " __local " << this->asString() << " locAlpha[" << localWorkgroupSize << "];" << std::endl;
        sourceStream << std::endl;
    }

    sourceStream << " " << this->asString() << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl
            << std::endl;
    sourceStream << " " << this->asString() << " myResult = 0.0;" << std::endl << std::endl;


    //caching data in register array, this also requires loading the data into the registers (in contrast using pointers to data directly)
    if (parameters["KERNEL_STORE_DATA"].compare("array") == 0) {
        for (size_t i = 0; i < dataBlockSize; i++) {
            sourceStream << indent << this->asString() << " data_" << i << "[" << dims << "];" << std::endl;
        }
        sourceStream << std::endl;
        for (size_t i = 0; i < dataBlockSize; i++) {
            for (size_t d = 0; d < dims; d++) {
                sourceStream << indent << getData(d, i) << " = ptrData[" << i << " + (" << dataBlockSize
                        << " * globalIdx) + (resultSize * " << d << ")];" << std::endl;
            }
            sourceStream << std::endl;
        }
    } else if (parameters["KERNEL_STORE_DATA"].compare("register") == 0) {
        for (size_t i = 0; i < dataBlockSize; i++) {
            for (size_t d = 0; d < dims; d++) {
                sourceStream << indent << this->asString() << " " << getData(d, i) << " = ptrData[" << i << " + ("
                        << dataBlockSize << " * globalIdx) + (resultSize * " << d << ")];" << std::endl;
            }
            sourceStream << std::endl;
        }
    }

    sourceStream << indent << "size_t dimLevelIndex;" << std::endl;

    sourceStream << std::endl;

    if (useLocalMemory) {
        sourceStream << " // Iterate over all grid points (fast ones, with cache)" << std::endl;
        sourceStream << " uint chunkSizeGrid = end_grid - start_grid;" << std::endl;
        sourceStream << " uint fastChunkSizeGrid = (chunkSizeGrid / " << localWorkgroupSize << ") * "
                << localWorkgroupSize << ";" << std::endl;
        sourceStream << " for(int j = start_grid; j < start_grid + fastChunkSizeGrid; j+=" << localWorkgroupSize << ")"
                << std::endl;
        sourceStream << " {" << std::endl;

        if (dims > maxDimUnroll) {
            sourceStream << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
            sourceStream << "   locLevel[(localIdx*" << dims << ")+ d] = ptrLevel[((j+localIdx)*" << dims << ")+ d];"
                    << std::endl;
            sourceStream << "   locIndex[(localIdx*" << dims << ")+ d] = ptrIndex[((j+localIdx)*" << dims << ")+ d];"
                    << std::endl;
            sourceStream << "}" << std::endl;
        } else {
            for (size_t d = 0; d < dims; d++) {
                sourceStream << "   locLevel[(localIdx*" << dims << ")+" << d << "] = ptrLevel[((j+localIdx)*" << dims
                        << ")+" << d << "];" << std::endl;
                sourceStream << "   locIndex[(localIdx*" << dims << ")+" << d << "] = ptrIndex[((j+localIdx)*" << dims
                        << ")+" << d << "];" << std::endl;
            }
        }

        sourceStream << "   locAlpha[localIdx] = ptrAlpha[j+localIdx];" << std::endl;
        sourceStream << "   barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
        sourceStream << std::endl;
        sourceStream << "   for(int k = 0; k < " << localWorkgroupSize << "; k++)" << std::endl;
        sourceStream << "   {" << std::endl;
        sourceStream << "     curSupport = locAlpha[k];" << std::endl << std::endl;
    } else {
        sourceStream << "   for(int k = start_grid; k < end_grid; k++)" << std::endl;
        sourceStream << "   {" << std::endl;
        sourceStream << "     curSupport = ptrAlpha[k];" << std::endl << std::endl;
    }

//    size_t evenDimensions = (dims / dataBlockSize) * dataBlockSize;
//    size_t remainingDimensions = evenDimensions - remainingDimensions;
//    sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, dims, "");

    if (dims > maxDimUnroll) {
        sourceStream << indent2 << "for (size_t unrollDim = 0; unrollDim < " << ((dims / maxDimUnroll) * maxDimUnroll)
                << "; unrollDim += " << maxDimUnroll << ") {" << std::endl;
        sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, std::min(maxDimUnroll, dims), "unrollDim");
        sourceStream << indent2 << "}" << std::endl;

        if (dims % maxDimUnroll != 0) {
            sourceStream
                    << this->unrolledBasisFunctionEvalulation(dims, (dims / maxDimUnroll) * maxDimUnroll, dims, "");
        }
    } else {
        sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, dims, "");
    }

    sourceStream << "     myResult += curSupport;" << std::endl;
    sourceStream << "  }" << std::endl;
    sourceStream << std::endl;

    if (useLocalMemory) {
        sourceStream << "  barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
        sourceStream << " }" << std::endl;
    }

    sourceStream << std::endl;
    sourceStream << " ptrResult[globalIdx] = myResult;" << std::endl;
    sourceStream << "}" << std::endl;

    if (parameters.getAsBoolean("WRITE_SOURCE")) {
        this->writeSource("StreamingOCLMP_mult.cl", sourceStream.str());
    }

    return sourceStream.str();
}

}
}

