/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLFastMultiPlattform/StreamingModOCLFastMultiPlatformKernelSourceBuilder.hpp>

//#include <StreamingOCLParameters.hpp>

namespace SGPP {
namespace datadriven {

std::string
StreamingModOCLFastMultiPlatformKernelSourceBuilder::generateSourceMult() {

  if (firstDeviceConfig["REUSE_SOURCE"].getBool()) {
    return this->reuseSource("streamingModOCLFastMP_mult.cl");
  }

  std::stringstream sourceStream;

  if (this->asString() == "double") {
    sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl <<
                 std::endl;
  }

  // write signature
  sourceStream << "__kernel" << std::endl;
  sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize <<
               ", 1, 1)))" << std::endl;
  sourceStream << "void multOCL(__global const " << this->asString() <<
               "* ptrLevel," << std::endl;
  sourceStream << "             __global const " << this->asString() <<
               "* ptrIndex," << std::endl;
  sourceStream << "             __global const " << this->asString() <<
               "* ptrData," << std::endl;
  sourceStream << "             __global const " << this->asString() <<
               "* ptrAlpha," << std::endl;
  sourceStream << "             __global       " << this->asString() <<
               "* ptrResult," << std::endl;
  sourceStream << "             uint resultSize," << std::endl;
  sourceStream << "             uint start_grid," << std::endl;
  sourceStream << "             uint end_grid) {" << std::endl;
  sourceStream << indent << "int globalIdx = get_global_id(0);" << std::endl;
  sourceStream << indent << "int localIdx = get_local_id(0);" << std::endl;
  sourceStream << std::endl;

  //blocked result variables
  for (size_t i = 0; i < dataBlockSize; i++) {
    sourceStream << indent << this->asString() << " curSupport_" << i << ";" <<
                 std::endl;
    sourceStream << indent << this->asString() << " myResult_" << i << " = 0.0;" <<
                 std::endl << std::endl;
  }

  //caching data in register array, this also requires loading the data into the registers (in contrast using pointers to data directly)
  if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("array") == 0) {
    for (size_t i = 0; i < dataBlockSize; i++) {
      sourceStream << indent << this->asString() << " data_" << i << "[" << dims <<
                   "];" << std::endl;
    }

    sourceStream << std::endl;

    for (size_t i = 0; i < dataBlockSize; i++) {
      for (size_t d = 0; d < dims; d++) {
        sourceStream << indent << getData(d,
                                          i) << " = ptrData[" << i << " + (" << dataBlockSize
                     << " * globalIdx) + (resultSize * " << d << ")];" << std::endl;
      }

      sourceStream << std::endl;
    }
  } else if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("register") ==
             0) {
    for (size_t i = 0; i < dataBlockSize; i++) {
      for (size_t d = 0; d < dims; d++) {
        sourceStream << indent << this->asString() << " " << getData(d,
                     i) << " = ptrData[" << i << " + ("
                     << dataBlockSize << " * globalIdx) + (resultSize * " << d << ")];" << std::endl;
      }

      sourceStream << std::endl;
    }
  }

  sourceStream << indent << "size_t dimLevelIndex;" << std::endl;

  sourceStream << std::endl;

  //arrays for prefetching a chunk of levels and indices in shared memory (in contrast to using pointer to level/index directly)
  if (useLocalMemory) {
    sourceStream << indent << "__local " << this->asString() << " locLevel[" <<
                 dims* localWorkgroupSize << "];"
                 << std::endl;
    sourceStream << indent << "__local " << this->asString() << " locIndex[" <<
                 dims* localWorkgroupSize << "];"
                 << std::endl;
    sourceStream << indent << "__local " << this->asString() << " locAlpha[" <<
                 localWorkgroupSize << "];"
                 << std::endl;
    sourceStream << std::endl;

    sourceStream << indent <<
                 "// grid is always divisible by workgroup size (by padding)" << std::endl;
    sourceStream << indent << "uint chunkSizeGrid = end_grid - start_grid;" <<
                 std::endl;
    sourceStream << indent << "uint fastChunkSizeGrid = (chunkSizeGrid / " <<
                 localWorkgroupSize << ") * "
                 << localWorkgroupSize << ";" << std::endl;
    sourceStream << indent <<
                 "for(int j = start_grid; j < start_grid + fastChunkSizeGrid; j+="
                 << localWorkgroupSize << ") {" << std::endl;

    for (size_t d = 0; d < dims; d++) {
      sourceStream << indent2 << "locLevel[(localIdx * " << dims << ") + " << d
                   << "] = ptrLevel[((j + localIdx) * " << dims << ") + " << d << "];" <<
                   std::endl;
      sourceStream << indent2 << "locIndex[(localIdx * " << dims << ") + " << d
                   << "] = ptrIndex[((j + localIdx) * " << dims << ") + " << d << "];" <<
                   std::endl;
    }

    sourceStream << indent2 << "locAlpha[localIdx] = ptrAlpha[j + localIdx];" <<
                 std::endl;
    sourceStream << indent2 << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
    sourceStream << std::endl;
    sourceStream << indent2 << "for (size_t m = 0; m < " << localWorkgroupSize <<
                 "; m++) {" << std::endl;

    for (size_t i = 0; i < dataBlockSize; i++) {
      sourceStream << indent3 << "curSupport_" << i << " = locAlpha[m];" << std::endl
                   << std::endl;
    }
  } else {
    sourceStream << indent <<
                 "// Iterate over all grid points (slow ones, without cache)" << std::endl;
    sourceStream << indent << "for (int m = start_grid; m < end_grid; m++) {" <<
                 std::endl;

    for (size_t i = 0; i < dataBlockSize; i++) {
      sourceStream << indent2 << "curSupport_" << i << " = ptrAlpha[m];" << std::endl
                   << std::endl;
    }
  }

  if (dims > maxDimUnroll) {
    sourceStream << indent2 << "for (size_t unrollDim = 0; unrollDim < " << ((
                   dims / maxDimUnroll) * maxDimUnroll)
                 << "; unrollDim += " << maxDimUnroll << ") {" << std::endl;
    sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0,
                 std::min(maxDimUnroll, dims), "unrollDim");

    sourceStream << indent2 << "}" << std::endl;

    if (dims % maxDimUnroll != 0) {
      sourceStream
          << this->unrolledBasisFunctionEvalulation(dims,
              (dims / maxDimUnroll) * maxDimUnroll, dims, "");
    }
  } else {
    sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, dims, "");
  }

  for (size_t i = 0; i < dataBlockSize; i++) {
    sourceStream << indent2 << "myResult_" << i << " += curSupport_" << i << ";" <<
                 std::endl;
  }

  if (useLocalMemory) {
    sourceStream << indent2 << "}" << std::endl;
    sourceStream << std::endl;
    sourceStream << indent2 << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
  }

  sourceStream << indent << "}" << std::endl;
  sourceStream << std::endl;

  for (size_t i = 0; i < dataBlockSize; i++) {
    sourceStream << indent << "ptrResult[(" << dataBlockSize << " * globalIdx) + "
                 << i << "] = myResult_" << i
                 << ";" << std::endl;
  }

  sourceStream << "}" << std::endl;

  if (firstDeviceConfig["WRITE_SOURCE"].getBool()) {
    this->writeSource("streamingModOCLFastMP_mult.cl", sourceStream.str());
  }

  return sourceStream.str();
}

}
}
