/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include <sstream>

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLMask/StreamingModOCLMaskKernelSourceBuilder.hpp>

namespace SGPP {
namespace datadriven {

std::string StreamingModOCLMaskKernelSourceBuilder::generateSourceMult() {

    if ((*parameters)["REUSE_SOURCE"].getBool()) {
        return this->reuseSource("streamingModOCLMask_mult.cl");
    }

    std::stringstream sourceStream;
    if (this->asString() == "double") {
        sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
    }

    sourceStream << "__kernel" << std::endl;
    sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))" << std::endl;
    sourceStream << "void multOCLMask(__global const " << this->asString() << "* ptrLevel," << std::endl;
    sourceStream << "           __global const " << this->asString() << "* ptrIndex," << std::endl;
    sourceStream << "           __global const " << this->asString() << "* ptrMask," << std::endl; // not needed for this kernel, but there for uniformity
    sourceStream << "           __global const " << this->asString() << "* ptrOffset," << std::endl; // not needed for this kernel, but there for uniformity
    sourceStream << "           __global const " << this->asString() << "* ptrData," << std::endl;
    sourceStream << "           __global const " << this->asString() << "* ptrAlpha," << std::endl;
    sourceStream << "           __global       " << this->asString() << "* ptrResult," << std::endl;
    sourceStream << "           uint resultSize," << std::endl;
    sourceStream << "           uint start_grid," << std::endl;
    sourceStream << "           uint end_grid) " << std::endl;
    sourceStream << "{" << std::endl;
    sourceStream << "   int globalIdx = get_global_id(0);" << std::endl;
    sourceStream << "   int localIdx = get_local_id(0);" << std::endl;
    sourceStream << std::endl;

    if (useLocalMemory) {
        sourceStream << "   __local " << this->asString() << " locLevel[" << dims * localWorkgroupSize << "];"
                << std::endl;
        sourceStream << "   __local " << this->asString() << " locIndex[" << dims * localWorkgroupSize << "];"
                << std::endl;
        sourceStream << "   __local " << this->asString() << " locMask[" << dims * localWorkgroupSize << "];"
                << std::endl;
        sourceStream << "   __local " << this->asString() << " locOffset[" << dims * localWorkgroupSize << "];"
                << std::endl;
        sourceStream << "   __local " << this->asString() << " locAlpha[" << localWorkgroupSize << "];" << std::endl;
        sourceStream << std::endl;
    }

    sourceStream << "   " << this->asString() << " eval, index_calc, abs, last, localSupport, curSupport;" << std::endl
            << std::endl;
    sourceStream << "   " << this->asString() << " myResult = ptrResult[globalIdx];" << std::endl << std::endl;
    sourceStream << "   // Create registers for the data" << std::endl;

    for (size_t d = 0; d < dims; d++) {
        sourceStream << " " << this->asString() << " data_" << d << " = ptrData[globalIdx+(resultSize*" << d << ")];"
                << std::endl;
    }

    sourceStream << std::endl;
    if (useLocalMemory) {
        sourceStream << "   // Iterate over all grid points (fast ones, with cache)" << std::endl;
        sourceStream << " uint chunkSizeGrid = end_grid - start_grid;" << std::endl;
        sourceStream << " uint fastChunkSizeGrid = (chunkSizeGrid / " << localWorkgroupSize << ") * "
                << localWorkgroupSize << ";" << std::endl;
        sourceStream << " for(int j = start_grid; j < start_grid + fastChunkSizeGrid; j+=" << localWorkgroupSize << ")"
                << std::endl;
        sourceStream << "   {" << std::endl;

        for (size_t d = 0; d < dims; d++) {
            sourceStream << "     locLevel[(localIdx*" << dims << ")+" << d << "] = ptrLevel[((j+localIdx)*" << dims
                    << ")+" << d << "];" << std::endl;
            sourceStream << "     locIndex[(localIdx*" << dims << ")+" << d << "] = ptrIndex[((j+localIdx)*" << dims
                    << ")+" << d << "];" << std::endl;
            sourceStream << "     locMask[(localIdx*" << dims << ")+" << d << "] = ptrMask[((j+localIdx)*" << dims
                    << ")+" << d << "];" << std::endl;
            sourceStream << "     locOffset[(localIdx*" << dims << ")+" << d << "] = ptrOffset[((j+localIdx)*" << dims
                    << ")+" << d << "];" << std::endl;
        }

        sourceStream << "       locAlpha[localIdx] = ptrAlpha[j+localIdx];" << std::endl;
        sourceStream << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
        sourceStream << std::endl;
        sourceStream << "       for(int k = 0; k < " << localWorkgroupSize << "; k++)" << std::endl;
        sourceStream << "       {" << std::endl;
        sourceStream << "           curSupport = locAlpha[k];" << std::endl << std::endl;

        for (size_t d = 0; d < dims; d++) {
            sourceStream << "         eval = ((locLevel[(k*" << dims << ")+" << d << "]) * (data_" << d << "));"
                    << std::endl;
            sourceStream << "         index_calc = eval - (locIndex[(k*" << dims << ")+" << d << "]);" << std::endl;
            sourceStream << "         abs = as_" << this->asString() << "(as_" << this->intAsString()
                    << "(index_calc) | as_" << this->intAsString() << "(locMask[(k*" << dims << ")+" << d << "]));"
                    << std::endl;
            sourceStream << "         last = locOffset[(k*" << dims << ")+" << d << "] + abs;" << std::endl;
            sourceStream << "         localSupport = fmax(last, 0.0" << this->constSuffix() << ");" << std::endl;
            sourceStream << "         curSupport *= localSupport;" << std::endl << std::endl;
        }

        sourceStream << "           myResult += curSupport;" << std::endl;
        sourceStream << "       }" << std::endl;
        sourceStream << std::endl;
        sourceStream << "       barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
        sourceStream << "   }" << std::endl;
        sourceStream << std::endl;
        sourceStream << "   // Iterate over all grid points (slow ones, without cache)" << std::endl;
        sourceStream << " for(int m = start_grid + fastChunkSizeGrid; m < end_grid; m++)" << std::endl;
        sourceStream << "   {" << std::endl;
        sourceStream << "       curSupport = ptrAlpha[m];" << std::endl << std::endl;

        for (size_t d = 0; d < dims; d++) {
            sourceStream << "     eval = ((ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << "));"
                    << std::endl;
            sourceStream << "     index_calc = eval - (ptrIndex[(m*" << dims << ")+" << d << "]);" << std::endl;
            sourceStream << "     abs = as_" << this->asString() << "(as_" << this->intAsString()
                    << "(index_calc) | as_" << this->intAsString() << "(ptrMask[(m*" << dims << ")+" << d << "]));"
                    << std::endl;
            sourceStream << "     last = ptrOffset[(m*" << dims << ")+" << d << "] + abs;" << std::endl;
            sourceStream << "     localSupport = fmax(last, 0.0" << this->constSuffix() << ");" << std::endl;
            sourceStream << "     curSupport *= localSupport;" << std::endl << std::endl;
        }

        sourceStream << "       myResult += curSupport;" << std::endl;
        sourceStream << "   }" << std::endl;
    } else {
        sourceStream << "   // Iterate over all grid points (without cache)" << std::endl;
        sourceStream << " for(int m = start_grid; m < end_grid; m++)" << std::endl;
        sourceStream << "   {" << std::endl;
        sourceStream << "       curSupport = ptrAlpha[m];" << std::endl << std::endl;

        for (size_t d = 0; d < dims; d++) {
            sourceStream << "     eval = ((ptrLevel[(m*" << dims << ")+" << d << "]) * (data_" << d << "));"
                    << std::endl;
            sourceStream << "     index_calc = eval - (ptrIndex[(m*" << dims << ")+" << d << "]);" << std::endl;
            sourceStream << "     abs = as_" << this->asString() << "(as_" << this->intAsString()
                    << "(index_calc) | as_" << this->intAsString() << "(ptrMask[(m*" << dims << ")+" << d << "]));"
                    << std::endl;
            sourceStream << "     last = ptrOffset[(m*" << dims << ")+" << d << "] + abs;" << std::endl;
            sourceStream << "     localSupport = fmax(last, 0.0" << this->constSuffix() << ");" << std::endl;
            sourceStream << "     curSupport *= localSupport;" << std::endl << std::endl;
        }

        sourceStream << "       myResult += curSupport;" << std::endl;
        sourceStream << "   }" << std::endl;

    }
    sourceStream << std::endl;
    sourceStream << "   ptrResult[globalIdx] = myResult;" << std::endl;
    sourceStream << "}" << std::endl;

    if ((*parameters)["WRITE_SOURCE"].getBool()) {
        this->writeSource("streamingModOCLMask_mult.cl", sourceStream.str());
    }

    return sourceStream.str();
}

}
}

