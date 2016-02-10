/*
 * OCLKernelBuilder.hpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/StreamingOCLMultiPlatformKernelSourceBuilder.hpp>

#include <fstream>
#include <sstream>

#include <sgpp/base/exception/operation_exception.hpp>

//#include <StreamingOCLMultiPlatformParameters.hpp>

namespace SGPP {
namespace datadriven {

std::string
StreamingOCLMultiPlatformKernelSourceBuilder::generateSourceMultTrans() {

  if (firstDeviceConfig["REUSE_SOURCE"].getBool()) {
    return this->reuseSource("StreamingOCLMultiPlatform_multTrans.cl");
  }

  size_t localWorkgroupSize = firstDeviceConfig["LOCAL_SIZE"].getUInt();
  bool useLocalMemory = firstDeviceConfig["KERNEL_USE_LOCAL_MEMORY"].getBool();
  uint64_t maxDimUnroll = firstDeviceConfig["KERNEL_MAX_DIM_UNROLL"].getUInt();

  std::stringstream sourceStream;

  if (this->asString() == "double") {
    sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl <<
                 std::endl;
  }

  sourceStream << "__kernel" << std::endl;

  sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize <<
               ", 1, 1)))" << std::endl;
  sourceStream << "void multTransOCL(__global const " << this->asString() <<
               "* ptrLevel," << std::endl;
  sourceStream << "                  __global const " << this->asString() <<
               "* ptrIndex," << std::endl;
  sourceStream << "                  __global const " << this->asString() <<
               "* ptrData," << std::endl;
  sourceStream << "                  __global const " << this->asString() <<
               "* ptrSource," << std::endl;
  sourceStream << "                  __global       " << this->asString() <<
               "* ptrResult," << std::endl;
  sourceStream << "                  uint sourceSize," << std::endl;
  sourceStream << "                  uint start_data," << std::endl;
  sourceStream << "                  uint end_data) {" << std::endl;
  sourceStream << indent << "int globalIdx = get_global_id(0);" << std::endl;
  sourceStream << indent << "int localIdx = get_local_id(0);" << std::endl;
  sourceStream << std::endl;

  for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
    sourceStream << indent << this->asString() << " myResult_" << gridIndex <<
                 " = 0.0;" << std::endl;
  }

  sourceStream << std::endl;

  if (useLocalMemory) {
    sourceStream << indent << "__local " << this->asString() << " locData[" << dims*
                 localWorkgroupSize << "];"
                 << std::endl;
    sourceStream << indent << "__local " << this->asString() << " locSource[" <<
                 localWorkgroupSize << "];"
                 << std::endl << std::endl;
  }

  // create a register storage for the level and index of the grid points of the work item
  if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("array") == 0) {
    for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
      sourceStream << indent << this->asString() << " level_" << gridIndex << "[" <<
                   dims << "];" << std::endl;
      sourceStream << indent << this->asString() << " index_" << gridIndex << "[" <<
                   dims << "];" << std::endl;

      for (size_t d = 0; d < dims; d++) {
        sourceStream << indent << "level_" << gridIndex << "[" << d << "] = ptrLevel[(("
                     << transGridBlockSize << " * globalIdx + "
                     << gridIndex << ") * " << dims << ")+" << d << "];" << std::endl;
        sourceStream << indent << "index_" << gridIndex << "[" << d << "] = ptrIndex[(("
                     << transGridBlockSize << " * globalIdx + "
                     << gridIndex << ") * " << dims << ")+" << d << "];" << std::endl;

      }

      sourceStream << std::endl;
    }

    sourceStream << std::endl;
  } else if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("register") ==
             0) {
    for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {

      for (size_t d = 0; d < dims; d++) {
        sourceStream << indent << this->asString() << " level_" << gridIndex << "_" << d
                     << " = ptrLevel[((" << transGridBlockSize << " * globalIdx + "
                     << gridIndex << ") * " << dims << ")+" << d << "];" << std::endl;
        sourceStream << indent << this->asString() << " index_" << gridIndex << "_" << d
                     << " = ptrIndex[((" << transGridBlockSize << " * globalIdx + "
                     << gridIndex << ") * " << dims << ")+" << d << "];" << std::endl;
      }

      sourceStream << std::endl;
    }

    sourceStream << std::endl;
  }

  sourceStream << indent << "// Iterate over all data points" << std::endl;

  if (useLocalMemory) {
    sourceStream << indent << "for(int i = start_data; i < end_data; i+=" <<
                 localWorkgroupSize << ") {"
                 << std::endl;

    if (dims > maxDimUnroll) {
      sourceStream << indent2 << "for (size_t d = 0; d < " << dims << "; d++) {" <<
                   std::endl;
      sourceStream << indent3 << "locData[(d * " << localWorkgroupSize
                   << ")+(localIdx)] = ptrData[(d * sourceSize) + (localIdx + i)];" << std::endl;
      sourceStream << indent2 << "}" << std::endl;
    } else {
      for (size_t d = 0; d < dims; d++) {
        sourceStream << indent2 << "locData[(" << d << "*" << localWorkgroupSize <<
                     ")+(localIdx)] = ptrData[("
                     << d << "*sourceSize)+(localIdx+i)];" << std::endl;
      }
    }

    sourceStream << indent2 << "locSource[localIdx] = ptrSource[i+localIdx];" <<
                 std::endl;
    sourceStream << indent2 << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl <<
                 std::endl;
    sourceStream << indent2 << "for(int k = 0; k < " << localWorkgroupSize <<
                 "; k++) {" << std::endl;

    for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
      sourceStream << indent3 << this->asString() << " curSupport_" << gridIndex <<
                   " = locSource[k];"
                   << std::endl;
    }
  } else {
    sourceStream << indent << "for(int k = start_data; k < end_data; k++) {" <<
                 std::endl;

    for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
      sourceStream << indent2 << this->asString() << " curSupport_" << gridIndex <<
                   " = ptrSource[k];"
                   << std::endl;
    }
  }

  sourceStream << std::endl;

  if (dims > maxDimUnroll) {
    sourceStream << indent2 << "for (size_t unrollDim = 0; unrollDim < " << ((
                   dims / maxDimUnroll) * maxDimUnroll)
                 << "; unrollDim += " << maxDimUnroll << ") {" << std::endl;

    sourceStream << this->unrolledBasisFunctionEvalulationTrans(dims, 0,
                 std::min(maxDimUnroll, dims), "unrollDim");
    sourceStream << indent2 << "}" << std::endl;

    if (dims % maxDimUnroll != 0) {
      sourceStream
          << this->unrolledBasisFunctionEvalulationTrans(dims,
              (dims / maxDimUnroll) * maxDimUnroll, dims,
              "");
    }

  } else {
    sourceStream << this->unrolledBasisFunctionEvalulationTrans(dims, 0, dims, "");
  }

  sourceStream << std::endl;

  for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
    sourceStream << indent2 << "myResult_" << gridIndex << " += curSupport_" <<
                 gridIndex << ";" << std::endl;
  }

  sourceStream << indent << "}" << std::endl << std::endl;

  if (useLocalMemory) {
    sourceStream << indent2 << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
    sourceStream << indent << "}" << std::endl;
  }

  for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {

    sourceStream << indent << "ptrResult[(" << transGridBlockSize <<
                 " * globalIdx) + " << gridIndex
                 << "] = myResult_" << gridIndex << ";" << std::endl;
  }

  sourceStream << "}" << std::endl;

  if (firstDeviceConfig["WRITE_SOURCE"].getBool()) {
    this->writeSource("StreamingOCLMultiPlatform_multTrans.cl", sourceStream.str());
  }

  return sourceStream.str();
}

}
}

