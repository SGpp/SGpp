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
StreamingModOCLFastMultiPlatformKernelSourceBuilder::generateSourceMultTrans() {

  if (firstDeviceConfig["REUSE_SOURCE"].getBool()) {
    return this->reuseSource("streamingModOCLFastMP_multTrans.cl");
  }

  size_t localWorkgroupSize = firstDeviceConfig["LOCAL_SIZE"].getUInt();
  size_t transGridBlockSize =
    firstDeviceConfig["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt();
  size_t transDataBlockSize =
    firstDeviceConfig["KERNEL_TRANS_DATA_BLOCK_SIZE"].getUInt();

  std::stringstream sourceStream;

  // preamble and signature
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
  sourceStream << "                  uint resultOffset," << std::endl;
  sourceStream << "                  uint sourceSize," << std::endl;
  sourceStream << "                  uint start_data," << std::endl;
  sourceStream << "                  uint end_data) {" << std::endl;
  //    stream_program_src << "   int globalSize = get_global_size(0);" << std::endl;
  sourceStream << indent << "int globalIdx = get_global_id(0);" << std::endl;
  sourceStream << indent << "int groupIdx = get_group_id(0);" << std::endl;
  sourceStream << indent << "int localIdx = get_local_id(0);" << std::endl;
  sourceStream << indent << "int assignedComponentBase = resultOffset + groupIdx;"
               << std::endl;
  sourceStream << std::endl;

  // array for local reduction
  sourceStream << indent << "__local " << this->asString() << " resultsTemp[" <<
               localWorkgroupSize << "];" << std::endl;
  sourceStream << std::endl;

  // blocked result variables
  for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
    sourceStream << indent << this->asString() << " myResult_" << gridIndex <<
                 " = 0.0;" << std::endl << std::endl;
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
                     << transGridBlockSize
                     << " * assignedComponentBase + " << gridIndex << ") * " << dims << ") +" << d <<
                     "];" << std::endl;
        sourceStream << indent << "index_" << gridIndex << "[" << d << "] = ptrIndex[(("
                     << transGridBlockSize
                     << " * assignedComponentBase + " << gridIndex << ") * " << dims << ") +" << d <<
                     "];" << std::endl;
      }

      sourceStream << std::endl;
    }

    sourceStream << std::endl;
  } else if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("register") ==
             0) {
    for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {


      for (size_t d = 0; d < dims; d++) {
        sourceStream << indent << this->asString() << " level_" << gridIndex << "_" << d
                     << " = ptrLevel[((" << transGridBlockSize
                     << " * assignedComponentBase + " << gridIndex << ") * " << dims << ") +" << d <<
                     "];" << std::endl;
        sourceStream << indent << this->asString() << " index_" << gridIndex << "_" << d
                     << " = ptrIndex[((" << transGridBlockSize
                     << " * assignedComponentBase + " << gridIndex << ") * " << dims << ") +" << d <<
                     "];" << std::endl;
      }

      sourceStream << std::endl;
    }

    sourceStream << std::endl;
  }

  sourceStream << indent << "size_t dimDataIndex;" << std::endl;

  if (firstDeviceConfig["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
    sourceStream << indent << "size_t dimLevelIndex;" << std::endl;
  }

  // iterate the data set and evaluate the basis functions
  sourceStream << indent <<
               "for(int k = start_data + localIdx; k < end_data; k += "
               << transDataBlockSize* localWorkgroupSize << ") {" << std::endl;

  for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
    for (size_t dataIndex = 0; dataIndex < transDataBlockSize; dataIndex++) {
      sourceStream << indent2 << this->asString() << " curSupport_" << gridIndex <<
                   "_" << dataIndex
                   << " = ptrSource[k + " << (localWorkgroupSize * dataIndex) << "];" << std::endl;
    }
  }

  sourceStream << std::endl;

  for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {

    if (dims > maxDimUnroll) {
      sourceStream << indent2 << "for (size_t unrollDim = 0; unrollDim < "
                   << ((dims / maxDimUnroll) * maxDimUnroll) << "; unrollDim += " << maxDimUnroll
                   << ") {"
                   << std::endl;
      sourceStream
          << this->unrolledBasisFunctionEvalulationTrans(dims, 0, std::min(maxDimUnroll,
              dims), "unrollDim",
              gridIndex);

      sourceStream << indent2 << "}" << std::endl;

      if (dims % maxDimUnroll != 0) {
        sourceStream
            << this->unrolledBasisFunctionEvalulationTrans(dims,
                (dims / maxDimUnroll) * maxDimUnroll, dims,
                "", gridIndex);
      }
    } else {
      sourceStream << this->unrolledBasisFunctionEvalulationTrans(dims, 0, dims, "",
                   gridIndex);
    }
  }

  for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
    for (size_t dataIndex = 0; dataIndex < transDataBlockSize; dataIndex++) {
      sourceStream << indent2 << "myResult_" << gridIndex << " += curSupport_" <<
                   gridIndex << "_" << dataIndex
                   << ";" << std::endl;
    }
  }

  sourceStream << indent << "}" << std::endl << std::endl;

  for (size_t gridIndex = 0; gridIndex < transGridBlockSize; gridIndex++) {
    if (gridIndex > 0) {
      sourceStream << indent << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl <<
                   std::endl;
    }

    sourceStream << indent << "resultsTemp[localIdx] = myResult_" << gridIndex <<
                 ";" << std::endl << std::endl;
    sourceStream << indent << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl <<
                 std::endl;

    sourceStream << indent << "if (localIdx == 0) {" << std::endl;
    sourceStream << indent2 << this->asString() << " overallResult = 0.0;" <<
                 std::endl;
    sourceStream << indent2 << "for (int i = 0; i < " << localWorkgroupSize <<
                 "; i++) {" << std::endl;
    sourceStream << indent3 << "overallResult += resultsTemp[i];" << std::endl;
    sourceStream << indent2 << "}" << std::endl;
    sourceStream << indent2 << "ptrResult[" << transGridBlockSize <<
                 " * assignedComponentBase + " << gridIndex
                 << "] = overallResult;" << std::endl;
    sourceStream << indent << "}" << std::endl;
  }

  sourceStream << "}" << std::endl;

  if (firstDeviceConfig["WRITE_SOURCE"].getBool()) {
    this->writeSource("streamingModOCLFastMP_multTrans.cl", sourceStream.str());
  }

  return sourceStream.str();
}

}
}
