/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include "StreamingOCLKernelSourceBuilder.hpp"

#include <fstream>
#include <sstream>

#include <sgpp/base/exception/operation_exception.hpp>

//#include "StreamingOCLParameters.hpp"

namespace SGPP {
namespace datadriven {

std::string StreamingOCLKernelSourceBuilder::generateSourceMultTrans(size_t dims) {

    if (parameters.getAsBoolean("REUSE_SOURCE")) {
        return this->reuseSource("StreamingOCL_multTrans.cl");
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
    sourceStream << "void multTransOCL(__global const " << this->asString() << "* ptrLevel," << std::endl;
    sourceStream << "           __global const " << this->asString() << "* ptrIndex," << std::endl;
    sourceStream << "           __global const " << this->asString() << "* ptrData," << std::endl;
    sourceStream << "           __global const " << this->asString() << "* ptrSource," << std::endl;
    sourceStream << "           __global       " << this->asString() << "* ptrResult," << std::endl;
    sourceStream << "           uint sourceSize," << std::endl;
    sourceStream << "           uint start_data," << std::endl;
    sourceStream << "           uint end_data)" << std::endl;
    sourceStream << "{" << std::endl;
    sourceStream << " int globalIdx = get_global_id(0);" << std::endl;
    sourceStream << " int localIdx = get_local_id(0);" << std::endl;
    sourceStream << std::endl;
    sourceStream << " " << this->asString() << " eval, index_calc, abs, last, localSupport, curSupport;"
            << std::endl << std::endl;
    sourceStream << " " << this->asString() << " myResult = ptrResult[globalIdx];" << std::endl << std::endl;

    if (useLocalMemory) {
        sourceStream << " __local " << this->asString() << " locData[" << dims * localWorkgroupSize << "];"
                << std::endl;
        sourceStream << " __local " << this->asString() << " locSource[" << localWorkgroupSize << "];"
                << std::endl << std::endl;
    }

    if (dims > maxDimUnroll) {
        sourceStream << " " << this->asString() << " level[" << dims << "];" << std::endl;
        sourceStream << " " << this->asString() << " index[" << dims << "];" << std::endl;
        sourceStream << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
        sourceStream << "  level[d] = ptrLevel[(globalIdx*" << dims << ") + d];" << std::endl;
        sourceStream << "  index[d] = ptrIndex[(globalIdx*" << dims << ") + d];" << std::endl;
        sourceStream << "}" << std::endl;
    } else {
        for (size_t d = 0; d < dims; d++) {
            sourceStream << " " << this->asString() << " level_" << d << " = ptrLevel[(globalIdx*" << dims << ")+"
                    << d << "];" << std::endl;
            sourceStream << " " << this->asString() << " index_" << d << " = ptrIndex[(globalIdx*" << dims << ")+"
                    << d << "];" << std::endl;
        }
    }

    sourceStream << std::endl;
    sourceStream << " // Iterate over all grid points" << std::endl;

    if (useLocalMemory) {
        sourceStream << " for(int i = start_data; i < end_data; i+=" << localWorkgroupSize << ")" << std::endl;
        sourceStream << " {" << std::endl;

        if (dims > maxDimUnroll) {
            sourceStream << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;
            sourceStream << "   locData[(d * " << localWorkgroupSize
                    << ")+(localIdx)] = ptrData[(d * sourceSize) + (localIdx + i)];" << std::endl;
            sourceStream << "}" << std::endl;
        } else {
            for (size_t d = 0; d < dims; d++) {
                sourceStream << "   locData[(" << d << "*" << localWorkgroupSize << ")+(localIdx)] = ptrData[("
                        << d << "*sourceSize)+(localIdx+i)];" << std::endl;
            }
        }

        sourceStream << "   locSource[localIdx] = ptrSource[i+localIdx];" << std::endl;
        sourceStream << "   barrier(CLK_LOCAL_MEM_FENCE);" << std::endl << std::endl;
        sourceStream << "   for(int k = 0; k < " << localWorkgroupSize << "; k++)" << std::endl;
        sourceStream << "   {" << std::endl;

        sourceStream << "     curSupport = locSource[k];" << std::endl << std::endl;
    } else {
        sourceStream << "   for(int k = start_data; k < end_data; k++)" << std::endl;
        sourceStream << "   {" << std::endl;
        sourceStream << "     curSupport = ptrSource[k];" << std::endl << std::endl;
    }

    if (dims > maxDimUnroll) {
        sourceStream << "for (size_t d = 0; d < " << dims << "; d++) {" << std::endl;

        if (useLocalMemory) {
            sourceStream << "     eval = (level[d] * (locData[(d * " << localWorkgroupSize << ")+k]));"
                    << std::endl;
        } else {
            sourceStream << "     eval = (level[d] * (ptrData[(d * sourceSize) + k]));" << std::endl;
        }

        sourceStream << "     index_calc = eval - index[d];" << std::endl;
        sourceStream << "     abs = fabs(index_calc);" << std::endl;
        sourceStream << "     last = 1.0" << this->constSuffix() << " - abs;" << std::endl;
        sourceStream << "     localSupport = fmax(last, 0.0" << this->constSuffix() << ");" << std::endl;
        sourceStream << "     curSupport *= localSupport;" << std::endl;
        sourceStream << "}" << std::endl;
    } else {
        for (size_t d = 0; d < dims; d++) {
            if (useLocalMemory) {
                sourceStream << "     eval = ((level_" << d << ") * (locData[(" << d << "*" << localWorkgroupSize
                        << ")+k]));" << std::endl;
            } else {
                sourceStream << "     eval = ((level_" << d << ") * (ptrData[(" << d << "*sourceSize)+k]));"
                        << std::endl;
            }

            sourceStream << "     index_calc = eval - (index_" << d << ");" << std::endl;
            sourceStream << "     abs = fabs(index_calc);" << std::endl;
            sourceStream << "     last = 1.0" << this->constSuffix() << " - abs;" << std::endl;
            sourceStream << "     localSupport = fmax(last, 0.0" << this->constSuffix() << ");" << std::endl;
            sourceStream << "     curSupport *= localSupport;" << std::endl;
        }
    }

    sourceStream << std::endl << "     myResult += curSupport;" << std::endl;
    sourceStream << "   }" << std::endl << std::endl;

    if (useLocalMemory) {
        sourceStream << "   barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
        sourceStream << " }" << std::endl;
    }

    sourceStream << " ptrResult[globalIdx] = myResult;" << std::endl;
    sourceStream << "}" << std::endl;

    if (parameters.getAsBoolean("WRITE_SOURCE")) {
        this->writeSource("StreamingOCL_multTrans.cl", sourceStream.str());
    }

    return sourceStream.str();
}

}
}

