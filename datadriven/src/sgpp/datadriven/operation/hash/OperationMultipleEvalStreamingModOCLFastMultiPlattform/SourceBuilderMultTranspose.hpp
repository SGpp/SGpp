// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <fstream>
#include <string>
#include <algorithm>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>

namespace sgpp {
namespace datadriven {
namespace StreamingModOCLFastMultiPlatform {

template <typename real_type>
class SourceBuilderMultTranspose : public base::KernelSourceBuilderBase<real_type> {
 private:
  std::shared_ptr<base::OCLDevice> device;

  json::Node &kernelConfiguration;

  size_t dims;

  size_t localWorkgroupSize;
  bool useLocalMemory;
  uint64_t maxDimUnroll;
  size_t transGridBlockSize;
  size_t transDataBlockSize;

  std::string getDataTrans(std::string dim, size_t dataBlockingIndex) {
    std::stringstream output;
    output << "ptrData[dimDataIndex + " << (localWorkgroupSize * dataBlockingIndex) << "]";
    return output.str();
  }

  std::string getLevelTrans(std::string dim, size_t gridBlockingIndex) {
    std::stringstream output;
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
      output << "level_" << gridBlockingIndex << "[" << dim << "]";
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
      output << "level_" << gridBlockingIndex << "_" << dim;
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
      output << "ptrLevel[dimLevelIndex]";
    } else {
      throw base::operation_exception(
          "OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
  }

  std::string getIndexTrans(std::string dim, size_t gridBlockingIndex) {
    std::stringstream output;
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
      output << "index_" << gridBlockingIndex << "[" << dim << "]";
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
      output << "index_" << gridBlockingIndex << "_" << dim;
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
      output << "ptrIndex[dimLevelIndex]";
    } else {
      throw base::operation_exception(
          "OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
    }
    return output.str();
  }

  std::string unrolledBasisFunctionEvalulationTrans(size_t dims, size_t startDim, size_t endDim,
                                                    std::string unrollVariable,
                                                    size_t gridBlockIndex) {
    size_t transDataBlockSize = kernelConfiguration["KERNEL_TRANS_DATA_BLOCK_SIZE"].getUInt();
    size_t transGridBlockSize = kernelConfiguration["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt();

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

      output << this->indent[2] << "dimDataIndex = "
             << "(" << pointerAccess << " * sourceSize) + k;" << std::endl;

      if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
        output << this->indent[2] << "dimLevelIndex = "
               << "((" << transGridBlockSize << " * groupIdx + " << gridBlockIndex << ") * " << dims
               << ") +" << pointerAccess << ";" << std::endl;
      }

      output << this->indent[2] << "if (" << getLevelTrans(dString, gridBlockIndex) << " == 2.0"
             << this->constSuffix() << ") {" << std::endl;

      output << this->indent[2] << "} else if (" << getIndexTrans(dString, gridBlockIndex)
             << " == 1.0" << this->constSuffix() << ") {" << std::endl;

      for (size_t dataIndex = 0; dataIndex < transDataBlockSize; dataIndex++) {
        output << this->indent[3] << "curSupport_" << gridBlockIndex << "_" << dataIndex
               << " *= max(2.0" << this->constSuffix() << " - ("
               << getLevelTrans(dString, gridBlockIndex) << " * "
               << getDataTrans(dString, dataIndex) << "), 0.0" << this->constSuffix() << ") ;"
               << std::endl;
      }

      output << this->indent[2] << "} else if (" << getIndexTrans(dString, gridBlockIndex)
             << " == (" << getLevelTrans(dString, gridBlockIndex) << " - 1.0" << this->constSuffix()
             << ")) {" << std::endl;

      for (size_t dataIndex = 0; dataIndex < transDataBlockSize; dataIndex++) {
        output << this->indent[3] << "curSupport_" << gridBlockIndex << "_" << dataIndex
               << " *= max((" << getLevelTrans(dString, gridBlockIndex) << " * "
               << getDataTrans(dString, dataIndex) << ") - "
               << getIndexTrans(dString, gridBlockIndex) << " + 1.0, 0.0);" << std::endl;
      }

      output << this->indent[2] << "} else {" << std::endl;

      for (size_t dataIndex = 0; dataIndex < transDataBlockSize; dataIndex++) {
        output << this->indent[3] << "curSupport_" << gridBlockIndex << "_" << dataIndex
               << " *= max(1.0" << this->constSuffix() << " - fabs("
               << getLevelTrans(dString, gridBlockIndex) << " * "
               << getDataTrans(dString, dataIndex) << " - "
               << getIndexTrans(dString, gridBlockIndex) << "), 0.0" << this->constSuffix() << ");"
               << std::endl;
      }

      output << this->indent[2] << "}" << std::endl;
    }
    return output.str();
  }

 public:
  SourceBuilderMultTranspose(std::shared_ptr<base::OCLDevice> device,
                             json::Node &kernelConfiguration, size_t dims)
      : device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
    localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
    localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    transGridBlockSize = kernelConfiguration["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt();
    transDataBlockSize = kernelConfiguration["KERNEL_TRANS_DATA_BLOCK_SIZE"].getUInt();
    maxDimUnroll = kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt();
  }

  std::string generateSource() {
    if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
      return this->reuseSource("StreamingModOCLFastMultiPlatform_multTrans.cl");
    }

    std::stringstream sourceStream;

    // preamble and signature
    //    if (this->floatType() == "double") {
    if (std::is_same<real_type, double>::value) {
      sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl
                   << std::endl;
    }

    sourceStream << "__kernel" << std::endl;
    sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))"
                 << std::endl;
    sourceStream << "void multTransOCL(__global const " << this->floatType() << "* ptrLevel,"
                 << std::endl;
    sourceStream << "                  __global const " << this->floatType() << "* ptrIndex,"
                 << std::endl;
    sourceStream << "                  __global const " << this->floatType() << "* ptrData,"
                 << std::endl;
    sourceStream << "                  __global const " << this->floatType() << "* ptrSource,"
                 << std::endl;
    sourceStream << "                  __global       " << this->floatType() << "* ptrResult,"
                 << std::endl;
    //    sourceStream << "                  uint resultOffset," << std::endl;
    sourceStream << "                  int sourceSize," << std::endl;
    sourceStream << "                  int start_data," << std::endl;
    sourceStream << "                  int end_data) {" << std::endl;
    //    stream_program_src << "   int globalSize = get_global_size(0);" << std::endl;
    sourceStream << this->indent[0] << "int globalIdx = get_global_id(0);" << std::endl;
    sourceStream << this->indent[0] << "int groupIdx = get_group_id(0);" << std::endl;
    sourceStream << this->indent[0] << "int localIdx = get_local_id(0);" << std::endl;
    sourceStream << this->indent[0] << "int assignedComponentBase = groupIdx;" << std::endl;
    //    sourceStream << this->indent[0] << "int assignedComponentBase = resultOffset + groupIdx;"
    //                 << std::endl;
    sourceStream << std::endl;

    // array for local reduction
    sourceStream << this->indent[0] << "__local " << this->floatType() << " resultsTemp["
                 << localWorkgroupSize << "];" << std::endl;
    sourceStream << std::endl;

    // blocked result variables
    for (size_t gridPoint = 0; gridPoint < transGridBlockSize; gridPoint++) {
      sourceStream << this->indent[0] << this->floatType() << " myResult_" << gridPoint << " = 0.0;"
                   << std::endl
                   << std::endl;
    }

    // create a register storage for the level and index of the grid points of the work item
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
      for (size_t gridPoint = 0; gridPoint < transGridBlockSize; gridPoint++) {
        sourceStream << this->indent[0] << this->floatType() << " level_" << gridPoint << "["
                     << dims << "];" << std::endl;
        sourceStream << this->indent[0] << this->floatType() << " index_" << gridPoint << "["
                     << dims << "];" << std::endl;

        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[0] << "level_" << gridPoint << "[" << d << "] = ptrLevel[(("
                       << transGridBlockSize << " * assignedComponentBase + " << gridPoint << ") * "
                       << dims << ") +" << d << "];" << std::endl;
          sourceStream << this->indent[0] << "index_" << gridPoint << "[" << d << "] = ptrIndex[(("
                       << transGridBlockSize << " * assignedComponentBase + " << gridPoint << ") * "
                       << dims << ") +" << d << "];" << std::endl;
        }
        sourceStream << std::endl;
      }
      sourceStream << std::endl;
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
      for (size_t gridPoint = 0; gridPoint < transGridBlockSize; gridPoint++) {
        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[0] << this->floatType() << " level_" << gridPoint << "_" << d
                       << " = ptrLevel[((" << transGridBlockSize << " * assignedComponentBase + "
                       << gridPoint << ") * " << dims << ") +" << d << "];" << std::endl;
          sourceStream << this->indent[0] << this->floatType() << " index_" << gridPoint << "_" << d
                       << " = ptrIndex[((" << transGridBlockSize << " * assignedComponentBase + "
                       << gridPoint << ") * " << dims << ") +" << d << "];" << std::endl;
        }
        sourceStream << std::endl;
      }
      sourceStream << std::endl;
    }

    sourceStream << this->indent[0] << "int dimDataIndex;" << std::endl;

    // TODO(pfandedd): what is this for?
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
      sourceStream << this->indent[0] << "int dimLevelIndex;" << std::endl;
    }

    // iterate the data set and evaluate the basis functions
    sourceStream << this->indent[0] << "for(int k = start_data + localIdx; k < end_data; k += "
                 << transDataBlockSize *localWorkgroupSize << ") {" << std::endl;

    for (size_t gridPoint = 0; gridPoint < transGridBlockSize; gridPoint++) {
      for (size_t dataIndex = 0; dataIndex < transDataBlockSize; dataIndex++) {
        sourceStream << this->indent[1] << this->floatType() << " curSupport_" << gridPoint << "_"
                     << dataIndex << " = ptrSource[k + " << (localWorkgroupSize * dataIndex) << "];"
                     << std::endl;
      }
    }
    sourceStream << std::endl;

    for (size_t gridPoint = 0; gridPoint < transGridBlockSize; gridPoint++) {
      if (dims > maxDimUnroll) {
        sourceStream << this->indent[1] << "for (int unrollDim = 0; unrollDim < "
                     << ((dims / maxDimUnroll) * maxDimUnroll) << "; unrollDim += " << maxDimUnroll
                     << ") {" << std::endl;
        sourceStream << this->unrolledBasisFunctionEvalulationTrans(
            dims, 0, std::min(maxDimUnroll, dims), "unrollDim", gridPoint);

        sourceStream << this->indent[1] << "}" << std::endl;
        if (dims % maxDimUnroll != 0) {
          sourceStream << this->unrolledBasisFunctionEvalulationTrans(
              dims, (dims / maxDimUnroll) * maxDimUnroll, dims, "", gridPoint);
        }
      } else {
        sourceStream << this->unrolledBasisFunctionEvalulationTrans(dims, 0, dims, "", gridPoint);
      }
    }

    for (size_t gridPoint = 0; gridPoint < transGridBlockSize; gridPoint++) {
      for (size_t dataIndex = 0; dataIndex < transDataBlockSize; dataIndex++) {
        sourceStream << this->indent[1] << "myResult_" << gridPoint << " += curSupport_"
                     << gridPoint << "_" << dataIndex << ";" << std::endl;
      }
    }

    sourceStream << this->indent[0] << "}" << std::endl
                 << std::endl;

    for (size_t gridPoint = 0; gridPoint < transGridBlockSize; gridPoint++) {
      if (gridPoint > 0) {
        sourceStream << this->indent[0] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                     << std::endl;
      }

      sourceStream << this->indent[0] << "resultsTemp[localIdx] = myResult_" << gridPoint << ";"
                   << std::endl
                   << std::endl;
      sourceStream << this->indent[0] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl
                   << std::endl;

      sourceStream << this->indent[0] << "if (localIdx == 0) {" << std::endl;
      sourceStream << this->indent[1] << this->floatType() << " overallResult = 0.0;" << std::endl;
      sourceStream << this->indent[1] << "for (int i = 0; i < " << localWorkgroupSize << "; i++) {"
                   << std::endl;
      sourceStream << this->indent[2] << "overallResult += resultsTemp[i];" << std::endl;
      sourceStream << this->indent[1] << "}" << std::endl;
      sourceStream << this->indent[1] << "ptrResult[" << transGridBlockSize
                   << " * assignedComponentBase + " << gridPoint << "] = overallResult;"
                   << std::endl;
      sourceStream << this->indent[0] << "}" << std::endl;
    }

    sourceStream << "}" << std::endl;

    if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
      this->writeSource("StreamingModOCLFastMultiPlatform_multTrans.cl", sourceStream.str());
    }

    return sourceStream.str();
  }
};
}  // namespace StreamingModOCLFastMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
