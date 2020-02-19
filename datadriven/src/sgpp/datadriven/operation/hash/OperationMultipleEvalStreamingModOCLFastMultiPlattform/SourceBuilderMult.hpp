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
class SourceBuilderMult : public base::KernelSourceBuilderBase<real_type> {
 private:
  std::shared_ptr<base::OCLDevice> device;

  json::Node &kernelConfiguration;

  size_t dims;

  size_t localWorkgroupSize;
  bool useLocalMemory;
  size_t dataBlockSize;
  size_t transGridBlockSize;
  uint64_t maxDimUnroll;
  size_t prefetchSize;

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
      throw base::operation_exception(
          "OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
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

      output << this->indent[2] << "dimLevelIndex = "
             << "(m * " << dims << ") + " << pointerAccess << ";" << std::endl;

      output << this->indent[2] << "if (" << levelAccess << " == 2.0" << this->constSuffix()
             << ") {" << std::endl;
      output << this->indent[2] << "} else if (" << indexAccess << " == 1.0" << this->constSuffix()
             << ") {" << std::endl;

      for (size_t i = 0; i < dataBlockSize; i++) {
        output << this->indent[3] << "curSupport_" << i << " *= max(2.0" << this->constSuffix()
               << " - (" << levelAccess << " * " << getData(dString, i) << "), 0.0"
               << this->constSuffix() << ") ;" << std::endl;
      }

      output << this->indent[2] << "} else if (" << indexAccess << " == " << levelAccess << " - 1.0"
             << this->constSuffix() << ")" << std::endl;
      output << this->indent[2] << "{" << std::endl;

      for (size_t i = 0; i < dataBlockSize; i++) {
        output << this->indent[3] << "curSupport_" << i << " *= max((" << levelAccess << " * "
               << getData(dString, i) << ") - " << indexAccess << " + 1.0" << this->constSuffix()
               << ", 0.0" << this->constSuffix() << ");" << std::endl;
      }

      output << this->indent[2] << "} else {" << std::endl;

      for (size_t i = 0; i < dataBlockSize; i++) {
        output << this->indent[3] << "curSupport_" << i << " *= max(1.0" << this->constSuffix()
               << " - fabs(" << levelAccess << " * " << getData(dString, i) << " - " << indexAccess
               << "), 0.0" << this->constSuffix() << ");" << std::endl;
      }
      output << this->indent[2] << "}" << std::endl;
    }
    return output.str();
  }

 public:
  SourceBuilderMult(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration,
                    size_t dims)
      : device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
    localWorkgroupSize = kernelConfiguration["LOCAL_SIZE"].getUInt();
    useLocalMemory = kernelConfiguration["KERNEL_USE_LOCAL_MEMORY"].getBool();
    dataBlockSize = kernelConfiguration["KERNEL_DATA_BLOCK_SIZE"].getUInt();
    transGridBlockSize = kernelConfiguration["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt();
    maxDimUnroll = kernelConfiguration["KERNEL_MAX_DIM_UNROLL"].getUInt();
    prefetchSize = kernelConfiguration["KERNEL_PREFETCH_SIZE"].getUInt();
  }

  std::string generateSource() {
    if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
      return this->reuseSource("StreamingModOCLFastMultiPlatform_mult.cl");
    }

    std::stringstream sourceStream;

    if (std::is_same<real_type, double>::value) {
      sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl
                   << std::endl;
    }

    // write signature
    sourceStream << "__kernel" << std::endl;
    sourceStream << "__attribute__((reqd_work_group_size(" << localWorkgroupSize << ", 1, 1)))"
                 << std::endl;
    sourceStream << "void multOCL(__global const " << this->floatType() << "* ptrLevel,"
                 << std::endl;
    sourceStream << "             __global const " << this->floatType() << "* ptrIndex,"
                 << std::endl;
    sourceStream << "             __global const " << this->floatType() << "* ptrData,"
                 << std::endl;
    sourceStream << "             __global const " << this->floatType() << "* ptrAlpha,"
                 << std::endl;
    sourceStream << "             __global       " << this->floatType() << "* ptrResult,"
                 << std::endl;
    sourceStream << "             int resultSize," << std::endl;
    sourceStream << "             int start_grid," << std::endl;
    sourceStream << "             int end_grid) {" << std::endl;
    sourceStream << this->indent[0] << "int globalIdx = get_global_id(0);" << std::endl;
    sourceStream << this->indent[0] << "int localIdx = get_local_id(0);" << std::endl;
    sourceStream << std::endl;

    // blocked result variables
    for (size_t i = 0; i < dataBlockSize; i++) {
      sourceStream << this->indent[0] << this->floatType() << " curSupport_" << i << ";"
                   << std::endl;
      sourceStream << this->indent[0] << this->floatType() << " myResult_" << i << " = 0.0;"
                   << std::endl
                   << std::endl;
    }

    // caching data in register array, this also requires loading the data into the registers (in
    // contrast using pointers to data directly)
    if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
      for (size_t i = 0; i < dataBlockSize; i++) {
        sourceStream << this->indent[0] << this->floatType() << " data_" << i << "[" << dims << "];"
                     << std::endl;
      }
      sourceStream << std::endl;
      for (size_t i = 0; i < dataBlockSize; i++) {
        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[0] << getData(d, i) << " = ptrData[" << i << " + ("
                       << dataBlockSize << " * globalIdx) + (resultSize * " << d << ")];"
                       << std::endl;
        }
        sourceStream << std::endl;
      }
    } else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
      for (size_t i = 0; i < dataBlockSize; i++) {
        for (size_t d = 0; d < dims; d++) {
          sourceStream << this->indent[0] << this->floatType() << " " << getData(d, i)
                       << " = ptrData[" << i << " + (" << dataBlockSize
                       << " * globalIdx) + (resultSize * " << d << ")];" << std::endl;
        }
        sourceStream << std::endl;
      }
    }

    sourceStream << this->indent[0] << "int dimLevelIndex;" << std::endl;

    sourceStream << std::endl;

    // arrays for prefetching a chunk of levels and indices in shared memory (in contrast to using
    // pointer to level/index directly)
    if (useLocalMemory) {
      sourceStream << this->indent[0] << "__local " << this->floatType() << " locLevel["
                   << dims * prefetchSize << "];" << std::endl;
      sourceStream << this->indent[0] << "__local " << this->floatType() << " locIndex["
                   << dims * prefetchSize << "];" << std::endl;
      sourceStream << this->indent[0] << "__local " << this->floatType() << " locAlpha["
                   << prefetchSize << "];" << std::endl;
      sourceStream << std::endl;

      sourceStream << this->indent[0] << "for(int j = start_grid; j < end_grid; j+=" << prefetchSize
                   << ") {" << std::endl;

      sourceStream << this->indent[1] << "if (localIdx < " << prefetchSize << ") {" << std::endl;

      for (size_t d = 0; d < dims; d++) {
        sourceStream << this->indent[1] << "locLevel[(localIdx * " << dims << ") + " << d
                     << "] = ptrLevel[((j + localIdx) * " << dims << ") + " << d << "];"
                     << std::endl;
        sourceStream << this->indent[1] << "locIndex[(localIdx * " << dims << ") + " << d
                     << "] = ptrIndex[((j + localIdx) * " << dims << ") + " << d << "];"
                     << std::endl;
      }

      sourceStream << this->indent[1] << "locAlpha[localIdx] = ptrAlpha[j + localIdx];"
                   << std::endl;

      sourceStream << this->indent[1] << "}" << std::endl;
      sourceStream << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
      sourceStream << std::endl;
      sourceStream << this->indent[1] << "for (int m = 0; m < " << prefetchSize << "; m++) {"
                   << std::endl;
      for (size_t i = 0; i < dataBlockSize; i++) {
        sourceStream << this->indent[2] << "curSupport_" << i << " = locAlpha[m];" << std::endl
                     << std::endl;
      }
    } else {
      sourceStream << this->indent[0] << "for (int m = start_grid; m < end_grid; m++) {"
                   << std::endl;
      for (size_t i = 0; i < dataBlockSize; i++) {
        sourceStream << this->indent[1] << "curSupport_" << i << " = ptrAlpha[m];" << std::endl
                     << std::endl;
      }
    }

    if (dims > maxDimUnroll) {
      sourceStream << this->indent[1] << "for (int unrollDim = 0; unrollDim < "
                   << ((dims / maxDimUnroll) * maxDimUnroll) << "; unrollDim += " << maxDimUnroll
                   << ") {" << std::endl;
      sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, std::min(maxDimUnroll, dims),
                                                             "unrollDim");

      sourceStream << this->indent[1] << "}" << std::endl;
      if (dims % maxDimUnroll != 0) {
        sourceStream << this->unrolledBasisFunctionEvalulation(
            dims, (dims / maxDimUnroll) * maxDimUnroll, dims, "");
      }
    } else {
      sourceStream << this->unrolledBasisFunctionEvalulation(dims, 0, dims, "");
    }

    for (size_t i = 0; i < dataBlockSize; i++) {
      sourceStream << this->indent[1] << "myResult_" << i << " += curSupport_" << i << ";"
                   << std::endl;
    }

    if (useLocalMemory) {
      sourceStream << this->indent[1] << "}" << std::endl;
      sourceStream << std::endl;
      sourceStream << this->indent[1] << "barrier(CLK_LOCAL_MEM_FENCE);" << std::endl;
    }
    sourceStream << this->indent[0] << "}" << std::endl;
    sourceStream << std::endl;

    for (size_t i = 0; i < dataBlockSize; i++) {
      sourceStream << this->indent[0] << "ptrResult[(" << dataBlockSize << " * globalIdx) + " << i
                   << "] = myResult_" << i << ";" << std::endl;
    }

    sourceStream << "}" << std::endl;

    if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
      this->writeSource("StreamingModOCLFastMultiPlatform_mult.cl", sourceStream.str());
    }

    return sourceStream.str();
  }
};
}  // namespace StreamingModOCLFastMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
